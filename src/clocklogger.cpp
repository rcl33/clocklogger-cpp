#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <cmath>
#include <csignal>

#include <sys/types.h>
#include <sys/stat.h>
#include <cerrno>
#include <cstring>

#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"

//#include "boost/filesystem.hpp"

using namespace std;
using namespace boost::gregorian;
using namespace boost::posix_time;
//namespace bfs = boost::filesystem;

#include "tickanalyst.h"
#include "daq.h"

//#define USE_TEST_DATA

sig_atomic_t recording = 1;

extern "C" void sigint_handler(int signum)
{
	assert (signum == SIGINT);
	recording = 0;
}

double symround(double x) {
	double y = std::floor( std::fabs(x) + 0.5 );
	return (x < 0.0)? -y : y;
}
	

void dumpSamples(vector<double> &y1, vector<double> &y2) {
	ofstream ofs("dumpSamples.txt", ifstream::out);
	size_t i, N;
	
	if (!ofs) return;
	
	N = y1.size();
	if (y2.size() < N) N = y2.size();
	
	//cout << "dumping " << N << " samples" << endl;
	
	for (i=0; i < N; i+=5) ofs << i << " " << y1[i] << " " << y2[i] << endl;
	
	//system("showDumpedSamples");
}

bool loadSamples(ifstream &ifs, int nSamples, vector<double> &PPS, vector<double> &IR) {
	double a, b;
	int i;
	PPS.clear(); IR.clear();
	
	if (!ifs) return false;
	
	for (i=0; (i < nSamples) && !ifs.eof(); i++) {
		ifs >> a >> b;
		IR.push_back(a-2.5); // swap order
		PPS.push_back(b-2.5);
	}
	return true;
}

void removeMeanValue(vector<double> &PPS, vector<double> &IR)
{
	double meanPPS=0, meanIR=0;
	vector<double>::iterator it;
	
	for (it=PPS.begin(); it!=PPS.end(); ++it) meanPPS += *it;
	for (it=IR.begin(); it!=IR.end(); ++it) meanIR += *it;
	
	meanPPS /= PPS.size();
	meanIR /= IR.size();
	
	for (it=PPS.begin(); it!=PPS.end(); ++it) *it -= meanPPS;
	for (it=IR.begin(); it!=IR.end(); ++it) *it -= meanIR;	
}

void openDataFile(ofstream &ofs, const date &day)
{
	mode_t perms = S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH; // permissions
	
	//date today = day_clock::local_day();
	ostringstream fname; fname.fill('0');
	
	// Check each level of data tree exists, create it if not.
	fname << "/home/clock/clocklogger/data";
	if ((mkdir(fname.str().c_str(),perms)==0) || (errno==EEXIST)) {
		
		fname << "/" << setw(4) << day.year();
		if ((mkdir(fname.str().c_str(),perms)==0) || (errno==EEXIST)) {

			fname << "/" << setw(2) << day.month().as_number();
			if ((mkdir(fname.str().c_str(),perms)==0) || (errno==EEXIST)) {
				
				fname << "/clock" << setw(4) << day.year() << "-"
					  << setw(2) << day.month().as_number() << "-"
				  	  << setw(2) << day.day().as_number() << ".txt";
				
				// Close file if it's open
				if (ofs.is_open()) ofs.close();
				
				// Open for appending, and configure to output in fixed point (no powers of 10)
				ofs.open(fname.str().c_str(), ofstream::app);
				ofs << fixed;
				
				cerr << "Opened new data file: " << fname.str() << endl;
				
				return;
			}
		}
	}
	
	cerr << "Couldn't mkdir " << fname.str() << ": " << strerror(errno) << endl;
}
	
void printUsage()
{
	cout<< "Usage: clocklogger [options]\n"
		<< "  -r         :  only outputs raw logging data to standard out\n"
		<< "  -d <drift> :  approximate current drift\n"
		<< "       (needed because pendulum only gives clock time to nearest second)\n"
		<< endl;
}

int main (int argc, char * const argv[]) {
	// Parse arguments
	int c;
	stringstream inarg;
	bool rawOutput = false;
	int secondsFast = 0;
	bool setSecsFast = false;
	double lastDrift = 0;
	while ((c = getopt (argc, argv, "hrd:")) != -1)
	{
		switch (c)
		{
			case 'h':
				printUsage();
				return 0;
			case 'r':
				rawOutput = true;
            	break;
			case 'd':
				inarg.str(optarg);
				inarg >> lastDrift;
				if (!inarg.fail()) // successfully got a double out of argument
					setSecsFast = true;
				else
					cerr << "Warning: Couldn't understand secondsFast number: '" << optarg << "'" << endl;
				break;
//			case '?':
//				cerr << "Unknown option '" << optopt << "'" << endl;
//             	return 1;
			default:
				abort ();
		}
	}
	
	// Actual program
	
	//cout << "ClockLogger v2.1" << endl;

	// Data acquisition
	vector<double> IR, PPS;                      // vectors to store data
	double sampleRate = 43478; //10000;          // sample rate
	TaskHandle taskHandle=0;                     // NIDAQmx task handle
	TickAnalyst analyst(PPS, IR, sampleRate);    // code to analyse data: gets samples from PPS & IR
	
	// Reserve space in PPS & IR to avoid reallocation (probably not necessary)
	PPS.reserve(6*sampleRate);
	IR.reserve(6*sampleRate);
	
	unsigned int pretrigger = 0.2 * sampleRate;	// target location of down tick in chunk
	unsigned int x0 = pretrigger; 			// position of down tick in samples (= pretrigger to get started)
	unsigned int nSamples;					// number of samples to get each time
	
	// Results
	double drift, amplitude, driftMod1;
	
	// Timing
	date currentDay(not_a_date_time);	// Today (initialise to invalid date to force data file opening)
	time_t recordtime;					// UNIX time when this tick happened
	ofstream datafile;					// Data file stream
	
	bool firstTime = true;				// flag for first time through loop
	
	
	// Set SIGINT handler to quit cleanly
	signal (SIGINT,sigint_handler);
	
	// Set up data input
#ifdef USE_TEST_DATA
	ifstream ifs("testData.txt", ifstream::in); // fake data...
#else	
	try {
		// Set up DAQ & start logging
		DAQSetUp(taskHandle, sampleRate);
		DAQStart(taskHandle);
	} catch (DAQError &e) {
		cout << "DAQ Error: " << e.what() << endl;
		DAQCleanUp(taskHandle);
		return -1;
	}
#endif
	
	if (!setSecsFast) {
		//cout << "How many seconds fast is the clock (rounded down)? ";
		cout << "What is the current drift? ";
		cin >> lastDrift;
		cout << endl;
	}
	
	if (!rawOutput) cout << " (you may see some signal errors while I'm locking onto the clock ticks...)" << endl;
	
	while (recording)
	{
		try
		{
			// Need to start a new data file?
			if ((true || !rawOutput) && (currentDay != day_clock::universal_day()))
			{
				// Close file if it's open
				//if (datafile.is_open()) datafile.close();
				
				// Open new data file
				currentDay = day_clock::universal_day();
				//dfname.str("");
				//dfname << "data2/" << currentDay.year() << "/"
				//	   << currentDay.month().as_number() << "/"
				//	   << currentDay.day().as_number();
				//cerr << "Trying to create: " << bfs::path(dfname.str()) << endl;
				//bfs::create_directory( bfs::path(dfname.str()) );

				//datafilename = "data/clockdata-" + to_iso_extended_string(currentDay) + ".txt";
				//datafile.open(datafilename.c_str(), ofstream::app);
				//datafile << fixed; // don't output in scientific format
				
				openDataFile (datafile, currentDay);
				
				if (!datafile) {
					cerr << "Could not open data file! \"" << datafile << "\"" << endl;
					DAQCleanUp(taskHandle);
					return -1;
				} else {
					//cerr << "Opened data file: " << datafilename << endl;
				}
			}
			
			// Decide how many samples to take to keep in sync with ticks
			if ( analyst.tickIndexKnown() ) x0 = analyst.tickIndex();
			nSamples = 3*sampleRate + ((signed int)x0 - (signed int)pretrigger)/3;
			//cout << "Recording 3 seconds " << showpos << (nSamples/sampleRate - 3) << noshowpos << endl;
			
			// Check this is reasonable
			if (nSamples > (3*sampleRate + 3*pretrigger)) {
				cerr << "Limiting recording length (too high)" << endl;
				nSamples = 3*sampleRate + 3*pretrigger;
			} else if (nSamples < (3*sampleRate - 3*pretrigger)) {
				cerr << "Limiting recording length (too low)" << endl;
				nSamples = 3*sampleRate - 3*pretrigger;
			}
			
#ifdef USE_TEST_DATA
			if (!loadSamples(ifs, nSamples, PPS, IR)) {
				cout << "out of data!" << endl; break;
			}
			usleep(500000);
#else
			DAQGetData(taskHandle, nSamples, PPS, IR); // get some data
#endif

			recordtime = ((int)time(NULL) / 3) * 3; // round down to ensure all times are 3*N seconds from midnight
			//recordtime = time(NULL); // -3; //+ (x0 - nSamples)/sampleRate;
			//removeMeanValue(PPS, IR);
			
			//cout << "got data: " << PPS.size() << " samples (" << PPS.size()/sampleRate << "s)" << endl;
			dumpSamples(PPS, IR); // dump samples to a file for debugging
			
			if (!rawOutput) cout << "[" << to_simple_string(from_time_t(recordtime)) << "] ";
			
			
			// Analyse this tick
			if ( !analyst.tick(driftMod1, amplitude) ) {
				cerr << "not locked onto ticks..." << endl;
				continue; // not locked
			}
			
			// Unwrap phase
			if (firstTime) {
				secondsFast = symround(-(driftMod1 - lastDrift)); // lastDrift has supplied value
				lastDrift = secondsFast + driftMod1; firstTime = false;
			}
			if ( (secondsFast+driftMod1 - lastDrift) > 0.5 ) {
				cerr << "(-1 sec: sF=" << secondsFast << ") ";
				--secondsFast;
			} else if ( (secondsFast+driftMod1 - lastDrift) < -0.5 ) {
				cerr << "(+1 sec: sF=" << secondsFast << ") ";
				++secondsFast;
			}
			lastDrift = drift = secondsFast + driftMod1;
			
			if (!rawOutput) cout << "drift = " << fixed << setprecision(6) << drift
								 << ", amp = " << fixed << setprecision(6) << amplitude << endl;
			
			// +/- 3 sec adjustment
			//drift += adjustClock;
			//secondsFast += adjustClock;
			//adjustClock = 0;
			
			// Write to file
			if (rawOutput)
				cout << fixed << setprecision(0) << recordtime << " "
					 << setprecision(6) << drift << " " << amplitude << endl;
			//else 
				datafile << setprecision(0) << recordtime << " "
						 << setprecision(6) << drift << " " << amplitude << endl;
			
		} catch (TickAnalyst::BadSignal &e) {
			// non-fatal error; carry on
			cerr << "Signal error: " << e.what();
			usleep(200000);
			cerr << " (...continuing...)" << endl;
		} catch (DAQError &e) {
			if (e.error_number == -200279) {
				cerr << "Warning: missed some samples... ";
			} else {
				cerr << "DAQ Error: " << e.what() << endl;
			}
			usleep(200000);
			cerr << "continuing..." << endl;
		} catch (runtime_error &e) {
			cerr << "Error: " << e.what() << endl;
			DAQCleanUp(taskHandle);
			return -1;
		}
	}
	
	cerr << "Cleaning up... ";
	DAQStop(taskHandle); // stop logging
	DAQCleanUp (taskHandle); // clean up (these do the same thing...)
	
	cerr << "done." << endl;
	
	return 0;
}

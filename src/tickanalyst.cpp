/*
 *  tickanalyst.cpp
 *  clocklogger
 *
 *  Created by Richard Lupton on 07/02/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <cmath>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "tickanalyst.h"

inline const double SQR(const double a) { return a*a; }

TickAnalyst::TickAnalyst(vector<double> &PPS, vector<double> &IR, double sampleRate)
: fs_(sampleRate)
, PPS_(PPS)
, IR_(IR)
, lastTickValid_(false)
, lastTickFromEnd_(4,0)
{
	startskip_  = int(fs_*0.060);   // ignore 60ms at start (arbitrary! from orig safetyMargin)
	pulsewidth_ = int(fs_*0.230);  // shim blocks IR for 230ms, ie time from -ve to +ve spike
	fitlength_  = int(fs_*0.030);
	wibbleskip_ = int(fs_*0.006);
	dcskip_     = int(fs_*0.008); // this must be enough to get before any mess at peak (eg. second peaks)
	dclength_   = int(fs_*0.015);
	
	ppsseek_    = int(fs_*0.100);  // look for PPS ±100ms from where we expect them
	
	ir2ndpulsewindow_ = int(fs_*2.2);
	
	irthreshold_  = 0.5; //0.25;
	ppsthreshold_ = 0.5; //0.15;
	intercept_    = 0.5;
	
	shimwidth_ = 24.0; // millimetres ~= milliradians (at 1m from pivot)
	
	downTick_ = 0;
	downTickKnowledge_ = 0;
	
	//PPS_ = NULL;
	//IR_ = NULL;
}

TickAnalyst::~TickAnalyst()
{
}


bool TickAnalyst::tick(double &drift, double &amplitude)
{
	vector<double> ppsx, irx;
	
	//PPS_ = PPS;
	//IR_ = IR;
	
	if (!findEdges(ppsx, irx)) {
		//cout << "Need to adjust synchronisation" << endl;
		return false;
	}
	
	analyseTick(ppsx, irx, drift, amplitude);
	return true;
}

/*void TickAnalyst::useChunk(const vector<double>& PPS, const vector<double>& IR)
{
	PPS_ = &PPS;
	IR_ = &IR;
}*/

// irx: 4 edges of 1 pendulum period
// ppsx: indices of as many PPS ticks as might be relevant
void TickAnalyst::analyseTick(const vector<double> &ppsx, const vector<double> &irx, double &drift, double &amplitude)
{
	int refppsx;
	unsigned int i;
	double mean_fs=0, var_fs=0, local_fs, period;
	vector<double> myppsx(ppsx);
	
	// Check we've got the right number of IR edges
	//if (irx.size() < 8) {
	//	std::ostringstream o; o << "Not enough IR edges found: got " << irx.size() << ", expected at least 8";
	//	throw ClockLogger::BadSignal(o.str());
	//}
	if (ppsx.size() < 2) {
		std::ostringstream o; o << "Not enough PPS edges found: got " << ppsx.size() << ", need at least 2";
		throw TickAnalyst::BadSignal(o.str());
	}	
	
	// Find mean sample rate from PPS signal
	for (i=0; i < (ppsx.size() - 1); i++) {
		mean_fs += (ppsx[i+1] - ppsx[i]);
		var_fs += SQR(ppsx[i+1] - ppsx[i]);
	}
	mean_fs /= (ppsx.size() - 1);
	var_fs = var_fs/(ppsx.size()-1) - SQR(mean_fs);
	if ( abs(mean_fs-fs_)/fs_ > 0.01 ) { // warn if sample rate error is above 1%
		std::ostringstream o; o << "Sample rate is off by too much: " << showpos << 1e6*abs(mean_fs-fs_)/fs_ << " ppm";
		throw TickAnalyst::BadSignal(o.str());
	}
	// PPS should be +/- 1us. Warn if std.deviation * 3, say, is greater than this.
	//  i.e. 9*variance > (1e-6 * fs)^2? but this is rather less than 1 sample...
	if ( var_fs > 200 ) { // warn if sample rate seems too variable (empirical)
		cerr << "Warning: PPS signal interval variance is high (" << int(var_fs) << ")." << endl;
		// Perhaps the signal is inverted and I'm measuring the falling rather than the rising edges?" << endl;
	}
	
	// Extrapolate the 1sec pulses at both ends in case the last IR is outside them
	//ppsx.insert( ppsx.begin(), *(ppsx.begin())-mean_fs );
	//ppsx.insert( ppsx.end(),   *(ppsx.last()) +mean_fs );
	//myppsx.push_back(ppsx.back() + mean_fs);
	
	// Find first PPS after down tick (irx[0])
	refppsx = -1;
	for (i=0; i < myppsx.size(); i++) {
		if ( myppsx[i] > irx[0] ) {
			refppsx = i; break;
		}
	}
	if ( refppsx < 0 ) {
		throw TickAnalyst::BadSignal("No PPS found after down tick!");
	}
		
	//local_fs = (refppsx > 0)? double(myppsx[refppsx] - myppsx[refppsx-1]) : mean_fs;
	if (refppsx > 0)
		local_fs = double(myppsx[refppsx] - myppsx[refppsx-1]);
	else if ((unsigned int)refppsx+1 < myppsx.size())
		local_fs = double(myppsx[refppsx+1] - myppsx[refppsx]);
	else
		local_fs = mean_fs;
	drift = (myppsx[refppsx] - irx[0]) / local_fs;
	
	if (drift > 1.5) { // must be a missing PPS, or something's gone wrong
		std::ostringstream o; o << setprecision(2) << fixed;
		o << "Time between down tick and next PPS too high: " << drift << "s";
		throw TickAnalyst::BadSignal(o.str());
	}
		
		
	// Period
	if (lastTickValid_) {
		for (period=0, i=0; i < 4; i++)
			period += irx[i] + lastTickFromEnd_[i];
		period /= (4*mean_fs);
	} else {
		period = 3; // best guess
	}
	
	for (i=0; i<4; i++)
		lastTickFromEnd_[i] = IR_.size() - irx[i];
	lastTickValid_ = true;
	
	// Amplitude
	amplitude = calcAmplitude(period, irx.begin(), mean_fs);
}

// Finds the edges in the ~tick's worth of data. Returns true if the chunk is
// aligned ok with the ticks (i.e. down tick first), false otherwise.
bool TickAnalyst::findEdges(dvec &ppsx, dvec &irx)
{
	unsigned int i;
	unsigned int x0, x0pps;
	vector<double> tempirx(4,-1);
	
	downTick_ = 0;
	downTickKnowledge_ = 0;
	
	// Find first pulse (after startskip)
	//  (if we can't find it, throw a NoEdgeInRange error)
	x0    = firstPositiveEdge (IR_,  startskip_, irthreshold_);
	x0pps = firstNegativeEdge (PPS_, startskip_, ppsthreshold_);
	
	downTick_ = x0; // might not be the correct (down) tick but better than nothing for tracking purposes
	++downTickKnowledge_;
	
	// Look for all 4 IR edges
	//  (if x0 is too late in data, ie x0+ir2ndpulsewindow would be out of
	//   bounds, firstXEdge will throw a NoEdgeInRange error. We catch this so
	//   we have a chance to return firstIR)
	
	if ( (x0+ir2ndpulsewindow_) >= IR_.size() ) {
		throw TickAnalyst::BadSignal("Not enough IR data after first tick");
	}
	
	try {
		//tempirx[0] = firstPositiveEdge ( IR_, x0-startskip_   , x0+1*pulsewidth_,     irthreshold_ );
		//tempirx[1] = firstNegativeEdge ( IR_, x0              , x0+2*pulsewidth_,     irthreshold_ );
		//tempirx[2] = firstPositiveEdge ( IR_, x0+2*pulsewidth_, x0+ir2ndpulsewindow_, irthreshold_ );
		//tempirx[3] = firstNegativeEdge ( IR_, x0+2*pulsewidth_, x0+ir2ndpulsewindow_, irthreshold_ );
		tempirx[0] = firstPositiveEdge ( IR_, x0-startskip_   , irthreshold_ );
		tempirx[1] = firstNegativeEdge ( IR_, x0              , irthreshold_ );
		tempirx[2] = firstPositiveEdge ( IR_, x0+2*pulsewidth_, irthreshold_ );
		tempirx[3] = firstNegativeEdge ( IR_, x0+2*pulsewidth_, irthreshold_ );
		//cerr << "\n>" << tempirx[0] << "\n>" << tempirx[1] << "\n>" << tempirx[2] << "\n>" << tempirx[3] << endl;
	} catch (NoEdgeInRange &e) {
		cerr << "findEdges: couldn't find all 4 IR edges";
		return false;
	}
	
	/* Find which is the 'down' tick: the IR signal looks like this:
	 *
	 * Tick:   ,u'     ,d'                    ,u'     ,d'
	 *
	 * (a):  | 0 1     2 3     (long gap)   |            
	 * (b):          | 0 1     (long gap)     2 3   |    
	 *
	 * The 3 seconds start either before an 'up' tick (a) or a 'down' tick (b).
	 * (down means pendulum is about to pass through centre, up is return swing)
	 * We distinguish the two cases by comparing the inner gap ([2] - [1]) with
	 * the outer gap ([0]+3secs - [3]).
	 */

	if ( (tempirx[2] - tempirx[1]) > ((tempirx[0] + (3*fs_)) -  tempirx[3]) ) {
		downTick_ = x0 = tempirx[0]; // 'down' tick is 1st
		++downTickKnowledge_;
	} else {
		downTick_ = x0 = tempirx[2]; // 'down' tick is 2nd
		return false; // need to align chunks with ticks before doing more
	}
		
	// and then precisely
	irx.clear();
	irx.push_back(findPrecisely(IR_, tempirx[0],  1));
	irx.push_back(findPrecisely(IR_, tempirx[1], -1));
	irx.push_back(findPrecisely(IR_, tempirx[2],  1));
	irx.push_back(findPrecisely(IR_, tempirx[3], -1));
		
	// PPS
	// how many 1s chunks fit after x0pps? (last need only be `ppsseek' long)
	//Npps = floor((double)(PPS_.size() - x0pps - ppsseek_)/fs_ + 1);
	//if (Npps < 1)
	//	throw ClockLogger::BadSignal("findPulses: No PPS signal found");
	
	ppsx.clear();  // already created by caller
	vector<double>::const_iterator ppsbeg = PPS_.begin();
	int x;
	
	for (i = 0; (x0pps+(i*fs_))+ppsseek_ < PPS_.size(); i++) {
		x = int(min_element(ppsbeg+(x0pps+(i*fs_))-ppsseek_, ppsbeg+(x0pps+(i*fs_))+ppsseek_) - ppsbeg);
		if ( PPS_.at(x) > -ppsthreshold_ )
			throw TickAnalyst::BadSignal("findPulses: Could not find PPS pulse");
		//ppsx[i] = x;
		ppsx.push_back(findPrecisely(PPS_, x, -1));
	}
	
	// Finally return index of first IR down tick
	//firstIR = x0;
	
	// Everything's ok if we got here
	return true;
}


/* Find 4 IR pulses, starting around x0: (w = pulsewidth)
 *  (1) +ve: x0-w    to   x0+w
 *  (2) -ve: x0      to   x0+2w
 *  (3) +ve: x0+2w   to   x0+ir2ndpulsewindow
 *  (4) -ve: x0+2w   to   x0+ir2ndpulsewindow
 */
/*bool TickAnalyst::find4ir(int x0, vector<double>& irx)
{
	vector<double>::const_iterator beg = IR_->begin();

	int x1 = int( max_element(beg + x0 - startskip_,    beg + x0 + 1*pulsewidth_)     - beg );
	int x2 = int( min_element(beg + x0             ,    beg + x0 + 2*pulsewidth_)     - beg );
	int x3 = int( max_element(beg + x0 + 2*pulsewidth_, beg + x0 + ir2ndpulsewindow_) - beg );
	int x4 = int( min_element(beg + x0 + 2*pulsewidth_, beg + x0 + ir2ndpulsewindow_) - beg );

	//cerr << x1 << " -> " << IR_->at(x1) << endl;
	//cerr << x2 << " -> " << IR_[x2] << endl;
	//cerr << x3 << " -> " << IR_[x3] << endl;
	//cerr << x4 << " -> " << IR_[x4] << endl;
	//cerr << IR_.size() << endl;
	
	if ( (IR_[x1] <  this->irthreshold_) ||
		 (IR_[x2] > -this->irthreshold_) ||
		 (IR_[x3] <  this->irthreshold_) ||
		 (IR_[x4] > -this->irthreshold_) ) {
		//cerr << x1 << " -> " << IR_[x1] << endl;
		//cerr << x2 << " -> " << IR_[x2] << endl;
		//cerr << x3 << " -> " << IR_[x3] << endl;
		//cerr << x4 << " -> " << IR_[x4] << endl;
		//return false;
		throw TickAnalyst::BadSignal("Find4IR: could not find 4 pulses");
	}
	
	irx.push_back(x1);
	irx.push_back(x2);
	irx.push_back(x3);
	irx.push_back(x4);
	
	return true;
}*/

double TickAnalyst::calcAmplitude(const double period, const vector<double>::const_iterator &it, const double fs)
{
	double om, num, den, t0;
	
	om = 2*M_PI/period / fs; // includes sample rate!
	
	num = sin( *(it+0)*om ) - sin( *(it+1)*om ) + sin( *(it+2)*om ) - sin( *(it+3)*om );
	den = cos( *(it+0)*om ) - cos( *(it+1)*om ) + cos( *(it+2)*om ) - cos( *(it+3)*om );
	t0 = atan2(num, den);
	
	return abs( shimwidth_ / ( sin(*(it+0)*om - t0) - sin(*(it+1)*om - t0) ) );
}

// Fits an exponential decay to the data, with peak roughly at x.
double TickAnalyst::findPrecisely(const vector<double> &y, const int &x, const int &sign) {
	double y0, a, b, siga, sigb, chi2;
	vector<double>::const_iterator it;
	
	vector<double> chunk( y.begin()+x+wibbleskip_,       y.begin()+x+fitlength_ );
	vector<double> dcoff( y.begin()+x-dcskip_-dclength_, y.begin()+x-dcskip_    );
	
	// zero before peak
	y0 = 0;
	for (it=dcoff.begin(); it != dcoff.end(); it++) 
		y0 += *it;
	y0 /= dcoff.size();
	
	int ysign = (sign > 0)? 1 : -1;
	
	//ofstream ofs ("data.txt", ifstream::out);
	//for (it=y.begin()+x-100; it != y.begin()+x+fitlength_; it++)
	//	ofs << *it << endl;

	
	fitExp(chunk, y0, ysign, a, b, siga, sigb, chi2);
	// returns y = a + bx (??)
	// so when y = threshold, x = (threshold - a) / b;
	return x+wibbleskip_ + (log(intercept_) - a) / b;
}

// Based on Numerical Recipies, p. 670.

void TickAnalyst::fitStraightLine(const dvec &x, const dvec &y, double &a, double &b, double &siga, double &sigb, double &chi2)
{
	int i;
	double t, sxoss, sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
	
	int ndata=x.size();
	b=0.0;
	
	// Accumulate sums
	for (i=0;i<ndata;i++) {
		sx += x[i];
		sy += y[i];
	}
	ss=ndata;
	sxoss=sx/ss;
	
	for (i=0;i<ndata;i++) {
		t=x[i]-sxoss;
		st2 += t*t;
		b += t*y[i];
	}
	
	// Solve for a, b, siga, sigb
	b /= st2;
	a=(sy-sx*b)/ss;
	siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
	sigb=sqrt(1.0/st2);
	
	// Calculate chi2 (for unweighted data evaluate typical sig using chi2, and adjust the standard deviations)
	chi2 = 0.0;
	//q=1.0;
	for (i=0; i<ndata;i++)
		chi2 += SQR(y[i]-a-b*x[i]);
	sigdat=sqrt(chi2/(ndata-2));
	siga *= sigdat;
	sigb *= sigdat;
}

// Modified version replaces y[i] with log(y[i]-y0) and assumes x = [0, 1, ... y.size()-1].
void TickAnalyst::fitExp(const dvec &y, const double &y0, const int &ysign, double &a, double &b, double &siga, double &sigb, double &chi2)
{
	int i;
	double t, sxoss, sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
	
	int ndata=y.size();
	b=0.0;
	

	
	// Accumulate sums
	for (i=0;i<ndata;i++) {
		sx += i;
		sy += log(ysign*(y[i]-y0));
	}
	ss=ndata;
	sxoss=sx/ss;
	
	for (i=0;i<ndata;i++) {
		t=i-sxoss;
		st2 += t*t;
		b += t*log(ysign*(y[i]-y0));
	}
	
	// Solve for a, b, siga, sigb
	b /= st2;
	a=(sy-sx*b)/ss;
	siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
	sigb=sqrt(1.0/st2);
	
	// Calculate chi2 (for unweighted data evaluate typical sig using chi2, and adjust the standard deviations)
	chi2 = 0.0;
	//q=1.0;
	for (i=0; i<ndata;i++)
		chi2 += SQR(log(ysign*(y[i]-y0))-a-b*i);
	sigdat=sqrt(chi2/(ndata-2));
	siga *= sigdat;
	sigb *= sigdat;
}

unsigned int TickAnalyst::firstPosNegEdge (const vector<double> &y, const unsigned int &start, const double &thresh, int sign) const
{
	vector<double>::const_iterator it(y.begin() + start);
	
	sign = (sign > 0)? 1 : -1;
	
	// check bounds
	if ( (it < y.begin()) ) throw TickAnalyst::NoEdgeInRange ();
	
	for ( ; it < y.end(); ++it ) {
		if ( (*it)*sign > thresh ) return int(it - y.begin());
		if ( isnan(*it) ) cerr << "Found NaN in raw input" << endl;
	}
	
	throw TickAnalyst::NoEdgeInRange();
}

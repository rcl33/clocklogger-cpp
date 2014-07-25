/*
 *  tickanalyst.h
 *  clocklogger
 *
 *  Created by Richard Lupton on 07/02/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <vector>
#include <algorithm>
#include <stdexcept>

using namespace std;

typedef vector<double> dvec;
typedef vector<int> ivec;

class TickAnalyst {
public:
	TickAnalyst(vector<double> &PPS, vector<double> &IR, double sampleRate);
	~TickAnalyst();
			
	// Analyse the recorded data in PPS & IR, returning drift & amplitude
	bool tick(double &drift, double &amplitude);

	// For synchronisation of recording to ticks: whether we know where a tick is, and if so, where
	bool tickIndexKnown() const { return downTickKnowledge_ > 0; }
	bool tickIsDownTick() const { return downTickKnowledge_ > 1; }
	unsigned int tickIndex() const { return downTick_; }

	class BadSignal : public std::runtime_error {
	public:
		BadSignal(const string &descr) : std::runtime_error(descr) { }
	};

	class NoEdgeInRange : public BadSignal {
	public:
		NoEdgeInRange() : BadSignal("No edge found in range") { }
	};
			
private:
	// Find the sample indices of the edges in the data
	bool findEdges(vector<double>& ppsx, vector<double>& irx);
	// Analyse these sample indices to find the drift and amplitude
	void analyseTick(const vector<double> &ppsx, const vector<double> &irx, double &drift, double &amplitude);

	//bool find4ir(int x0, vector<double>& irx);		
	
	// Calculate the amplitude of the pendulum
	double calcAmplitude(const double period, const vector<double>::const_iterator &it, const double fs);
			
	// Find the sub-sample location of an edge by fitting an exponential to the decay
	double findPrecisely(const vector<double> &y, const int &x, const int &sign);
	void fitStraightLine(const vector<double> &x, const vector<double> &y, double &a, double &b, double &siga, double &sigb, double &chi2);
	void fitExp(const vector<double> &y, const double &y0, const int &sign, double &a, double &b, double &siga, double &sigb, double &chi2);
	
	// Find the first +/- edge above the given threshold.	
	unsigned int firstPosNegEdge(const vector<double> &y, const unsigned int &start, const double &thresh, int sign) const;
	unsigned int firstPositiveEdge(const vector<double> &y, const unsigned int &start, const double &thresh) const { return firstPosNegEdge (y, start, thresh,  1); }
	unsigned int firstNegativeEdge(const vector<double> &y, const unsigned int &start, const double &thresh) const { return firstPosNegEdge (y, start, thresh, -1); }
	
	// Timing parameters
	unsigned int startskip_, pulsewidth_, fitlength_, wibbleskip_, dcskip_, dclength_, ir2ndpulsewindow_, ppsseek_;
	double intercept_, irthreshold_, ppsthreshold_, shimwidth_;	
			
	int fs_;
	const vector<double> &PPS_;
	const vector<double> &IR_;
			
	unsigned int downTick_;
	int downTickKnowledge_;
			
	bool lastTickValid_;
	vector<double> lastTickFromEnd_;
};

//~ inline
//~ unsigned int firstPositiveEdge (const vector<double>::iterator &start, const vector<double>::iterator &end, const double &thresh)
//~ {
	//~ vector<double>::const_iterator it(start);
	//~ for ( ; it < end; ++it )
		//~ if ( *it > thresh ) return int(it - start);
	//~ throw TickAnalyst::NoEdgeInRange();
//~ }

//~ inline
//~ unsigned int firstNegativeEdge (const vector<double>::const_iterator &start, const vector<double>::const_iterator &end, const double &thresh)
//~ {
	//~ vector<double>::const_iterator it(start);
	//~ for ( ; it < end; ++it )
		//~ if ( *it < -thresh ) return int(it - start);
	//~ throw TickAnalyst::NoEdgeInRange();
//~ }

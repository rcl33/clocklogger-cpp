/*
 *  testcases.cpp
 *  clocklogger
 *
 *  Created by Richard Lupton on 08/02/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

//#include "testcases.h"

#include <iostream>
#include <fstream>

using namespace std;

void testFind4IR (vector<double> &PPS, vector<double> &IR);
void testFindPulses (vector<double> &PPS, vector<double> &IR);

int main()
{
	double s1, s2;
	vector<double> PPS;
	vector<double> IR;
	PPS.reserve(397000);
	IR.reserve(397000);
	
	ifstream ifs ("test/chunk.txt", ifstream::in);
	if (!ifs) {
  		cout << "Cannot open test data file test/chunk.txt" << endl;
  		return -1;
  	}
	
	while (ifs >> s1 >> s2) {
		this->PPS.push_back(s1);
		this->IR.push_back(s2);
	}
}

void testFind4IR(vector<double> &PPS, vector<double> &IR)
{	
	ClockLogger cl;
	cl.useChunk(PPS, IR);
	
	int irAnswer[] = {103875, 113795, 178939, 188907, 236153, 246076, 311219, 321184};
	
	vector<double> irx(0);
	
	if (!cl.find4ir(irAnswer[0], irx)) {
		cerr << "problem in testFind4IR ()!" << endl;
		return;
	}
	
	int i;
	for (i = 0; i < 4; i++) {
		if ((int)irx[i] != irAnswer[i]) cerr << "IR #" << i << "doesn't match" << endl;
	}
}

void testcases::testFindPulses(vector<double> &PPS, vector<double> &IR)
{
	ClockLogger cl;
	cl.useChunk(PPS, IR);
	
	int ppsAnswer[] = {41332, 85428, 129524, 173620, 217715, 261811, 305907, 350003, 394098};
	int irAnswer[]  = {103875, 113795, 178939, 188907, 236153, 246076, 311219, 321184};
	
	vector<double> irx, ppsx;
	int x0;
	
	if (!cl.findPulses(ppsx, irx, x0)) {
		cerr << "problem in testFindPulses ()!" << endl;
		return;
	}

	
	int i;
	//cout << showpos;
	//cout << "PPS: ";
	if ( (irx.size() != 8) || (ppsx.size() != 8) ) {
		cout << "Wrong number of edges found: PPS " << ppsx.size()
			 << ", IR " << irx.size() << endl;
		return;
	}
	
	for ( i=0; i < irx.size(); i++ ) {
		if ((int)irx[i] != irAnswer[i]) {
			cerr << "Error: IR " << i+1 << " of 8 not found correctly: got " << irx[i] << ", expected " << irAnswer[i] << endl;
		}
	}
	for ( i=0; i < ppsx.size(); i++ ) {
		if ((int)ppsx[i] != ppsAnswer[i]) {
			cerr << "Error: pulse " << i+1 << " of 8 not found correctly: got " << ppsx[i] << ", expected " << ppsAnswer[i] << endl;
		}
	}
}

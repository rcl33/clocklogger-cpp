/*
 *  daq.h
 *  clocklogger
 *
 *  Created by Richard Lupton on 26/02/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <vector>
#include <stdexcept>

extern "C" {
#include <NIDAQmx.h>

int32 CVICALLBACK DoneCallback(TaskHandle taskHandle, int32 status, void *callbackData);
}
	
void DAQSetUp(TaskHandle &taskHandle,double sampleRate);
void DAQCleanUp(TaskHandle &taskHandle);
	
void DAQStart(TaskHandle &taskHandle);
void DAQStop(TaskHandle &taskHandle);
	
void DAQGetData(TaskHandle &taskHandle, const unsigned int &nSamples, std::vector<double> &chan1, std::vector<double> &chan2);
	
class DAQError : public std::runtime_error {
	public:
	DAQError(const std::string &descr, int32 err=0) : std::runtime_error(descr), error_number(err) { }
	int error_number;
};


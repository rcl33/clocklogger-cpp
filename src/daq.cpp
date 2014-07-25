/*********************************************************************
*
* ANSI C Example program:
*    ContAcq-IntClk.c
*
* Example Category:
*    AI
*
* Description:
*    This example demonstrates how to acquire a continuous amount of
*    data using the DAQ device's internal clock.
*
* Instructions for Running:
*    1. Select the physical channel to correspond to where your
*       signal is input on the DAQ device.
*    2. Enter the minimum and maximum voltage range.
*    Note: For better accuracy try to match the input range to the
*          expected voltage level of the measured signal.
*    3. Set the rate of the acquisition. Also set the Samples per
*       Channel control. This will determine how many samples are
*       read each time the while loop iterates. This also determines
*       how many points are plotted on the graph each iteration.
*    Note: The rate should be at least twice as fast as the maximum
*          frequency component of the signal being acquired.
*
* Steps:
*    1. Create a task.
*    2. Create an analog input voltage channel.
*    3. Set the rate for the sample clock. Additionally, define the
*       sample mode to be continuous.
*    4. Call the Start function to start the acquistion.
*    5. Read the data in a loop until the stop button is pressed or
*       an error occurs.
*    6. Call the Clear Task function to clear the task.
*    7. Display an error if any.
*
* I/O Connections Overview:
*    Make sure your signal input terminal matches the Physical
*    Channel I/O control. For further connection information, refer
*    to your hardware reference manual.
*
*********************************************************************/

#include <stdio.h>

#include <vector>

#include "daq.h"

#define DAQmxErrChk(functionCall) if( DAQmxFailed(error=(functionCall)) ) goto Error; else

//int32 CVICALLBACK EveryNCallback(TaskHandle taskHandle, int32 everyNsamplesEventType, uInt32 nSamples, void *callbackData);
int32 CVICALLBACK DoneCallback(TaskHandle taskHandle, int32 status, void *callbackData);

void DAQSetUp(TaskHandle &taskHandle, double sampleRate)
{
	int32       error=0;
	char        errBuff[2048]={'\0'};

	/*********************************************/
	// DAQmx Configure Code
	/*********************************************/
	DAQmxErrChk (DAQmxCreateTask("",&taskHandle));
	DAQmxErrChk (DAQmxCreateAIVoltageChan(taskHandle,"Dev1/ai0,Dev1/ai1","",DAQmx_Val_NRSE,-5.0,5.0,DAQmx_Val_Volts,NULL));
	DAQmxErrChk (DAQmxCfgSampClkTiming(taskHandle,"",sampleRate,DAQmx_Val_Rising,DAQmx_Val_ContSamps,9*sampleRate)); // 9s buffer -- bit on safe side

	//DAQmxErrChk (DAQmxRegisterEveryNSamplesEvent(taskHandle,DAQmx_Val_Acquired_Into_Buffer,1000,0,EveryNCallback,NULL));
	DAQmxErrChk (DAQmxRegisterDoneEvent(taskHandle,0,DoneCallback,NULL));
	return;
	
Error:
	if( DAQmxFailed(error) ) {
		DAQmxGetExtendedErrorInfo(errBuff,2048);
		throw DAQError(errBuff, error);
	}
}

void DAQCleanUp(TaskHandle &taskHandle)
{
	if( taskHandle != 0 ) {
		/*********************************************/
		// DAQmx Stop Code
		/*********************************************/
		DAQmxStopTask(taskHandle);
		DAQmxClearTask(taskHandle);
		taskHandle = 0;
	}
}

void DAQStart(TaskHandle &taskHandle)
{
	int32       error=0;
	char        errBuff[2048]={'\0'};
		
	if ( taskHandle != 0) {
		/*********************************************/
		// DAQmx Start Code
		/*********************************************/
		error = DAQmxStartTask(taskHandle);
	}
	
	if( DAQmxFailed(error) ) {
		DAQmxGetExtendedErrorInfo(errBuff,2048);
		throw DAQError(errBuff, error);
	}
}

void DAQStop(TaskHandle &taskHandle)
{
	//same as DAQCleanUp ()...
}

//int32 CVICALLBACK EveryNCallback(TaskHandle taskHandle, int32 everyNsamplesEventType, uInt32 nSamples, void *callbackData)
//{
void DAQGetData(TaskHandle &taskHandle, const unsigned int &nSamples, std::vector<double> &chan1, std::vector<double> &chan2)
{
	int32       error=0;
	char        errBuff[2048]={'\0'};
	float64     timeout;
	int32       read=0;
	
	if ( taskHandle == 0 ) return;
	
	float64 *data = new float64[2*nSamples];

	timeout = 9; // appropriate timeout?
	
	/*********************************************/
	// DAQmx Read Code
	/*********************************************/
	error = DAQmxReadAnalogF64(taskHandle,nSamples,timeout,DAQmx_Val_GroupByChannel,data,2*nSamples,&read,NULL);
	
	if( DAQmxFailed(error) ) {
		DAQmxGetExtendedErrorInfo(errBuff,2048);
		delete[] data;
		throw DAQError(errBuff, error);
	}	
	
	if( read>0 ) {
		chan1.assign( data, data+read );
		chan2.assign( data+read, data+2*read );
	}
	
	delete[] data;
}

int32 CVICALLBACK DoneCallback(TaskHandle taskHandle, int32 status, void *callbackData)
{
	char    errBuff[2048]={'\0'};

	// Check to see if an error stopped the task.
	if( DAQmxFailed(status) ) {
		DAQmxGetExtendedErrorInfo(errBuff,2048);
		DAQmxClearTask(taskHandle);
		printf("DAQmx Error: %s\n",errBuff);
	}
	return 0;
}

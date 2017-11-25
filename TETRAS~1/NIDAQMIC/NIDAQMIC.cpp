// NIDAQMIC.cpp: 콘솔 응용 프로그램의 진입점을 정의합니다.
//

#include "stdafx.h"
#include <NIDAQmx.h>

#define CHANNEL_PER_BUFFER 30000 // 절대 바꾸지 말 것
#define SAMPLING_RATE 25600.0
#define SAMPLE_PER_CHANNEL 25600
#define CHANNEL_COUNT 4

// NI-DAQ 장비 세팅
#define DAQmxErrChk(functionCall) if( DAQmxFailed(error=(functionCall)) ) goto Error; else

// NI 콜백
int32 CVICALLBACK EveryNCallback(TaskHandle taskHandle, int32 everyNsamplesEventType, uInt32 nSamples, void *callbackData);
int32 CVICALLBACK DoneCallback(TaskHandle taskHandle, int32 status, void *callbackData);
FILE* record;

int main(int argc, char* argv[])
{
	int32       error = 0;
	TaskHandle  taskHandle = 0;
	char        errBuff[2048] = { '\0' };

	// DAQmx Configure Code
	DAQmxErrChk(DAQmxCreateTask("", &taskHandle));
	DAQmxErrChk(DAQmxCreateAIMicrophoneChan(taskHandle, "cDAQ1Mod1/ai0", "", DAQmx_Val_PseudoDiff, DAQmx_Val_Pascals, 50, 120.0, DAQmx_Val_Internal, 0.004, NULL));
	DAQmxErrChk(DAQmxCfgSampClkTiming(taskHandle, "", SAMPLING_RATE, DAQmx_Val_Rising, DAQmx_Val_ContSamps, SAMPLE_PER_CHANNEL));
	DAQmxErrChk(DAQmxRegisterEveryNSamplesEvent(taskHandle, DAQmx_Val_Acquired_Into_Buffer, 100, 0, EveryNCallback, NULL));
	DAQmxErrChk(DAQmxRegisterDoneEvent(taskHandle, 0, DoneCallback, NULL));

	record = fopen("C:\\Users\\pblwo\\Desktop\\2017-10-14\\record.txt", "w");

	// DAQmx Start Code
	DAQmxErrChk(DAQmxStartTask(taskHandle));

	printf("Acquiring samples continuously. Press Enter to interrupt\n");
	getchar();

Error:
	if (DAQmxFailed(error))
		DAQmxGetExtendedErrorInfo(errBuff, 2048);

	if (taskHandle != 0)
	{
		// DAQmx Stop Code
		DAQmxStopTask(taskHandle);
		DAQmxClearTask(taskHandle);
	}

	if (DAQmxFailed(error))
		printf("DAQmx Error: %s\n", errBuff);

	fclose(record);
	printf("End of program, press Enter key to quit\n");
	getchar();
	return 0;
}

int32 CVICALLBACK EveryNCallback(TaskHandle taskHandle, int32 everyNsamplesEventType, uInt32 nSamples, void *callbackData)
{
	int32       i = 0;
	int32       error = 0;
	char        errBuff[2048] = { '\0' };
	static int  totalRead = 0;
	int32       read = 0;
	float64     data[CHANNEL_COUNT * CHANNEL_PER_BUFFER] = { 0.0, };  // 4채널 데이터 버퍼

	// DAQmx Read Code
	DAQmxErrChk(DAQmxReadAnalogF64(taskHandle, CHANNEL_PER_BUFFER, 10.0, DAQmx_Val_GroupByChannel,
		data, CHANNEL_COUNT * CHANNEL_PER_BUFFER, &read, NULL));

	if (read >= CHANNEL_PER_BUFFER)
	{
		printf("Acquired %d samples. Total %d\n", (int)read, (int)(totalRead += read));

		// 파일 저장
		for (i = 0; i < CHANNEL_PER_BUFFER; i++)
			fprintf(record, "%lf\t%lf\t%lf\t%lf\n", data[i], data[i + CHANNEL_PER_BUFFER], data[i + 2 * CHANNEL_PER_BUFFER], data[i + 3 * CHANNEL_PER_BUFFER]);

		fflush(stdout);
	}

Error:
	if (DAQmxFailed(error))
	{
		DAQmxGetExtendedErrorInfo(errBuff, 2048);

		// DAQmx Stop Code
		DAQmxStopTask(taskHandle);
		DAQmxClearTask(taskHandle);
		printf("DAQmx Error: %s\n", errBuff);
	}
	return 0;
}

int32 CVICALLBACK DoneCallback(TaskHandle taskHandle, int32 status, void *callbackData)
{
	int32   error = 0;
	char    errBuff[2048] = { '\0' };

	// Check to see if an error stopped the task.
	DAQmxErrChk(status);

Error:
	if (DAQmxFailed(error))
	{
		DAQmxGetExtendedErrorInfo(errBuff, 2048);
		DAQmxClearTask(taskHandle);
		printf("DAQmx Error: %s\n", errBuff);
	}
	return 0;
}

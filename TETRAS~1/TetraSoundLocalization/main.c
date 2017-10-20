/*
	�ۼ��� : ������
	�ۼ��� : 2017-09-14
	���� : �����ü ��ġ�� ����ũ���� ��̿� ���� ������ TDOA�� �̿��Ͽ� ������ ������ �� ������ �˾Ƴ��� �ڵ��
	AAI�˰����� �Ϻο� ������ �ڻ��� �˰����� ȥ�� �� �����Ͽ� �����Ͽ���.
	�Ʒ� �ڵ�� ���� �� 6���� TDOA�� �ִ� ���Ϸκ��� TDOA�� ��� ����Ѵ�.
*/

#pragma once

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <io.h>
#include <conio.h>
#include <string.h>

#include "sort.h"
#include "CPSP.h"
#include "tetra_sound_localization.h"

#define CHANNEL_PER_BUFFER 30000 // ���� �ٲ��� �� ��
#define SAMPLING_RATE 25600.0
#define SAMPLE_PER_CHANNEL 25600

//#define READ_FROM_DAQ
#define READ_FROM_FILE

FILE *record = NULL, *angle = NULL;

#ifdef READ_FROM_FILE

void main()
{
	double data[CHANNEL_COUNT * CHANNEL_PER_BUFFER] = { 0.0, };  // 4ä�� ������ ����
	double TDOA_01 = 0, TDOA_20 = 0, TDOA_30 = 0, TDOA_12 = 0, TDOA_13 = 0, TDOA_23 = 0;// Ÿ�� ������

	// OSH, AAI ���� ��
	double osh_azimuth[5] = { 0, }, osh_elevation[5] = { 0, };
	double aai_azimuth[5] = { 0, }, aai_elevation[5] = { 0, };
	double s_theta_all[10], s_piangle_all[10];

	// cpsp ����
	double MAX_Delay = 40, temp = 999;

	int i = 0, case_cnt = 1;

	record = fopen("C:\\Users\\pblwo\\Desktop\\����_3m_a000��_e30��_25600hz.txt", "r");
//	angle = fopen("./data/2017-10-09/noise/angle/awgn20_����_3m_a0��_e0��_25600hz.txt", "w");

	while (fscanf(record, "%lf\t%lf\t%lf\t%lf\n",
		&data[i], &data[i + CHANNEL_PER_BUFFER], &data[i + 2 * CHANNEL_PER_BUFFER], &data[i + 3 * CHANNEL_PER_BUFFER]) != EOF)
	{
		i++;
		if (i >= CHANNEL_PER_BUFFER)
		{
			fprintf(stdout, "%d\n", case_cnt);

			// 0. CPSP�� Ÿ�� ������ ���
			TDOA_01 = CPSP_FILT_normal(&data[0], &data[CHANNEL_PER_BUFFER], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;
			TDOA_12 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER], &data[CHANNEL_PER_BUFFER * 2], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;
			TDOA_20 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER * 2], &data[0], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;
			TDOA_30 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER * 3], &data[0], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;
			TDOA_13 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER], &data[CHANNEL_PER_BUFFER * 3], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;
			TDOA_23 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER * 2], &data[CHANNEL_PER_BUFFER * 3], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;

			do_osh_algorithm(osh_azimuth, osh_elevation, TDOA_01, TDOA_20, TDOA_30, TDOA_12, TDOA_13, TDOA_23);
			cord3_trans(osh_azimuth, osh_elevation, s_theta_all, s_piangle_all, 20, mic_distance);

			do_AAI_algorithm2(aai_azimuth, aai_elevation, TDOA_01, TDOA_20, TDOA_30, TDOA_12, TDOA_13, TDOA_23);
			cord3_trans(aai_azimuth, aai_elevation, s_theta_all, s_piangle_all, 20, mic_distance);

			//fprintf(stdout, "%lf\t%lf\t%lf\t%lf\n", osh_azimuth[4], osh_elevation[4], aai_azimuth[4], aai_elevation[4]);
			//fprintf(angle, "%lf\t%lf\t%lf\t%lf\n", osh_azimuth[4], osh_elevation[4], aai_azimuth[4], aai_elevation[4]);
			i = 0;
			case_cnt++;
			getchar();
		}
	}

	fclose(record);
//	fclose(angle);

}

#else READ_FROM_DAQ

#include <NIDAQmx.h>

// NI-DAQ ��� ����
#define DAQmxErrChk(functionCall) if( DAQmxFailed(error=(functionCall)) ) goto Error; else

// NI �ݹ�
int32 CVICALLBACK EveryNCallback(TaskHandle taskHandle, int32 everyNsamplesEventType, uInt32 nSamples, void *callbackData);
int32 CVICALLBACK DoneCallback(TaskHandle taskHandle, int32 status, void *callbackData);

int main(int argc, char* argv[])
{
	int32       error = 0;
	TaskHandle  taskHandle = 0;
	char        errBuff[2048] = { '\0' };

	// DAQmx Configure Code
	DAQmxErrChk(DAQmxCreateTask("", &taskHandle));
	DAQmxErrChk(DAQmxCreateAIMicrophoneChan(taskHandle, "cDAQ1Mod4/ai0:3", "", DAQmx_Val_PseudoDiff, DAQmx_Val_Pascals, 50, 120.0, DAQmx_Val_Internal, 0.004, NULL));
	DAQmxErrChk(DAQmxCfgSampClkTiming(taskHandle, "", SAMPLING_RATE, DAQmx_Val_Rising, DAQmx_Val_ContSamps, SAMPLE_PER_CHANNEL));
	DAQmxErrChk(DAQmxRegisterEveryNSamplesEvent(taskHandle, DAQmx_Val_Acquired_Into_Buffer, 100, 0, EveryNCallback, NULL));
	DAQmxErrChk(DAQmxRegisterDoneEvent(taskHandle, 0, DoneCallback, NULL));
		
	record = fopen("./data/2017-10-12/record/����_3m_a0��_e60��_25600hz.txt", "w");
	angle = fopen("./data/2017-10-12/angle/����_3m_a0��_e60��_25600hz.txt", "w");

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
	fclose(angle);

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
	float64     data[CHANNEL_COUNT * CHANNEL_PER_BUFFER] = { 0.0, };  // 4ä�� ������ ����

	// ���� ������, Ÿ�� ������
	double TDOA_01 = 0, TDOA_20 = 0, TDOA_30 = 0, TDOA_12 = 0, TDOA_13 = 0, TDOA_23 = 0;

	// OSH, AAI ���� ��
	double osh_azimuth[5] = { 0, }, osh_elevation[5] = { 0, };
	double aai_azimuth[5] = { 0, }, aai_elevation[5] = { 0, };
	double s_theta_all[10], s_piangle_all[10];

	// cpsp ����
	double MAX_Delay = 40, temp = 999;

	// �˰��� ���� ������ ���� ����
	static double case_cnt = 0, osh_error, aai_error = 0;
	static clock_t osh_start = 0, aai_start = 0, start = 0;
	static double osh_runtime = 0.0, aai_runtime = 0.0, runtime= 0.0;

	// DAQmx Read Code
	DAQmxErrChk(DAQmxReadAnalogF64(taskHandle, CHANNEL_PER_BUFFER, 10.0, DAQmx_Val_GroupByChannel, 
		data, CHANNEL_COUNT * CHANNEL_PER_BUFFER, &read, NULL));
	
	if (read >= CHANNEL_PER_BUFFER)
	{
		// printf("Acquired %d samples. Total %d\n", (int)read, (int)(totalRead += read));

		// ���� ����
		for (i = 0; i < CHANNEL_PER_BUFFER; i++)
			fprintf(record, "%lf\t%lf\t%lf\t%lf\n", data[i], data[i + CHANNEL_PER_BUFFER], data[i + 2 * CHANNEL_PER_BUFFER], data[i + 3 * CHANNEL_PER_BUFFER]);

		start = clock();

		// 0. CPSP�� Ÿ�� ������ ���
		TDOA_01 = CPSP_FILT_normal(&data[0], &data[CHANNEL_PER_BUFFER], gs_buff_cpsp, MAX_Delay) / SAMPLING_RATE;
		TDOA_12 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER], &data[CHANNEL_PER_BUFFER * 2], gs_buff_cpsp, MAX_Delay) / SAMPLING_RATE;
		TDOA_20 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER * 2], &data[0], gs_buff_cpsp, MAX_Delay) / SAMPLING_RATE;
		TDOA_30 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER * 3], &data[0], gs_buff_cpsp, MAX_Delay) / SAMPLING_RATE;
		TDOA_13 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER], &data[CHANNEL_PER_BUFFER * 3], gs_buff_cpsp, MAX_Delay) / SAMPLING_RATE;
		TDOA_23 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER * 2], &data[CHANNEL_PER_BUFFER * 3], gs_buff_cpsp, MAX_Delay) / SAMPLING_RATE;
		
		printf("total case : %2lf\n", case_cnt);

		// OSH
		//osh_start = clock();
		do_osh_algorithm(osh_azimuth, osh_elevation, TDOA_01, TDOA_20, TDOA_30, TDOA_12, TDOA_13, TDOA_23);
		cord3_trans(osh_azimuth, osh_elevation, s_theta_all, s_piangle_all, 20, mic_distance);		
		//osh_runtime += (double)(clock() - osh_start) / CLOCKS_PER_SEC; // ���� �ð� ����
	 
		//aai_start = clock();

		// AAI
		 do_AAI_algorithm2(aai_azimuth, aai_elevation, TDOA_01, TDOA_20, TDOA_30, TDOA_12, TDOA_13, TDOA_23);
		 cord3_trans(aai_azimuth, aai_elevation, s_theta_all, s_piangle_all, 20, mic_distance);

		//runtime += ((double)(clock() - start) / CLOCKS_PER_SEC);

		case_cnt++;

		fprintf(angle, "%lf\t%lf\t%lf\t%lf\n", osh_azimuth[4], osh_elevation[4], aai_azimuth[4], aai_elevation[4]);
		printf("%lf %lf %lf %lf\n", osh_azimuth[4], osh_elevation[4], aai_azimuth[4], aai_elevation[4]);
		//printf("%lf %lf\n", aai_azimuth[4], aai_elevation[4]);
		//printf("runtime : %lf\n", runtime);

		//printf("angle - osh : [%lf , %lf], aai : [%lf, %lf]\n", osh_azimuth[4], osh_elevation[4], aai_azimuth[4], aai_elevation[4]);
		//printf("runtime - osh : %lfms, aai : %lfms\n", osh_runtime, aai_runtime);
		//printf("average runtime - osh : %lfms, aai : %lfms\n", osh_runtime / case_cnt, aai_runtime / case_cnt);
		//printf("error rate - osh : %lf\%, aai : %lf\%\n", (osh_error / case_cnt) * 100, (aai_error / case_cnt) * 100);
		printf("=====================================\n\n");
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
#endif

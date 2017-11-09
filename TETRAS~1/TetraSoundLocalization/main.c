/*
	작성자 : 최지수
	작성일 : 2017-09-14
	설명 : 정사면체 배치의 마이크로폰 어레이에 들어온 음원의 TDOA를 이용하여 음원의 방위각 및 고도각을 알아내는 코드로
	AAI알고리즘의 일부와 오상헌 박사의 알고리즘을 혼용 및 배합하여 적용하였다.
	아래 코드는 라인 당 6개의 TDOA가 있는 파일로부터 TDOA를 얻어 사용한다.
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

#define CHANNEL_PER_BUFFER 30000 // 절대 바꾸지 말 것
#define SAMPLING_RATE 25600.0
#define SAMPLE_PER_CHANNEL 25600

//#define READ_FROM_DAQ
#define READ_FROM_FILE

FILE *record = NULL, *angle = NULL, *TDOA = NULL;

#ifdef READ_FROM_FILE

void main()
{
	srand(time(NULL));

	double data[CHANNEL_COUNT * CHANNEL_PER_BUFFER] = { 0.0, };  // 4채널 데이터 버퍼
	double TDOA_01 = 0, TDOA_20 = 0, TDOA_30 = 0, TDOA_12 = 0, TDOA_13 = 0, TDOA_23 = 0;// 타임 딜레이

	double azimuth[5] = { -999, -999, -999, -999, -999 }, elevation[5] = { -999, -999, -999, -999, -999 };
	// OSH, AAI 방위 고도
	double osh_azimuth[5] = { -999, -999, -999, -999, -999 }, osh_elevation[5] = { -999, -999, -999, -999, -999 };
	double aai_azimuth[5] = { -999, -999, -999, -999, -999 }, aai_elevation[5] = { -999, -999, -999, -999, -999 };
	double s_theta_all[10], s_piangle_all[10];

	// cpsp 인자
	double MAX_Delay = 40, temp = 999;
	double TDOA[6] = { 0.0 };
	int i = 0, case_cnt = 1;
	clock_t start = 0;
	double GT_azimuth = 0.0, GT_elevation = 0.0;
	double error_azimuth = 0.0, error_elevation = 0.0;
	double deg2rad = (M_PI / 180.0);
	double rad2deg = (180.0 / M_PI);
	record = fopen("C:\\Users\\pblwo\\Desktop\\noise\\고도 30도\\awgn20_녹음_6m_a060도_e30도_25600hz.txt", "r");
	if (record == NULL)
	{
		fprintf(stdout, "no file\n");
		return 0;
	}
	angle = fopen("C:\\Users\\pblwo\\Desktop\\noise\\고도 30도\\angle\\awgn20_각도_6m_a060도_e30도_25600hz.txt", "w");
	
//	TDOA = fopen("C:\\Users\\pblwo\\Desktop\\2017-10-14\\2017-10-14\\TDOA\\고도각 0도\\TDOA_3m_a000도_e00도_25600hz.txt", "w");
 
//while (fscanf(record, "%lf\t%lf\t%lf\t%lf\n",	
	//	&data[i], &data[i + CHANNEL_PER_BUFFER], &data[i + 2 * CHANNEL_PER_BUFFER], &data[i + 3 * CHANNEL_PER_BUFFER]) != EOF)
	
	while(case_cnt <= 1000)
	{
		rand_TDOA(TDOA, mic_distance, &GT_azimuth, &GT_elevation);
		TDOA_01 = TDOA[0];
		TDOA_12 = TDOA[1];
		TDOA_20 = TDOA[2];
		TDOA_30 = TDOA[3];
		TDOA_13 = TDOA[4];
		TDOA_23 = TDOA[5];
		
		//do_osh_algorithm(azimuth, elevation, TDOA_01, TDOA_20, TDOA_30, TDOA_12, TDOA_13, TDOA_23);		
		do_AAI_algorithm2(azimuth, elevation, TDOA_01, TDOA_20, TDOA_30, TDOA_12, TDOA_13, TDOA_23);
		cord3_trans(azimuth, elevation, s_theta_all, s_piangle_all, 20, mic_distance);
		//printf("%d : %lf %lf\n", case_cnt, GT_azimuth * rad2deg, GT_elevation * rad2deg);
		//printf("%d : %lf %lf\n", case_cnt, azimuth[4], elevation[4]);
		
		if (isnan(azimuth[4]) || isnan(elevation[4]))
			continue;

		azimuth[4] *= deg2rad;
		elevation[4] *= deg2rad;
		
		error_azimuth = fabs(((GT_azimuth - azimuth[4]) * rad2deg));
		error_elevation = fabs(((GT_elevation - elevation[4]) * rad2deg));
		case_cnt++;

		//if (i >= CHANNEL_PER_BUFFER)
		//{
		//	fprintf(stdout, "%d\n", case_cnt);

		//	// 0. CPSP로 타임 딜레이 계산
		//	TDOA_01 = CPSP_FILT_normal(&data[0], &data[CHANNEL_PER_BUFFER], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;
		//	TDOA_12 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER], &data[CHANNEL_PER_BUFFER * 2], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;
		//	TDOA_20 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER * 2], &data[0], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;
		//	TDOA_30 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER * 3], &data[0], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;
		//	TDOA_13 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER], &data[CHANNEL_PER_BUFFER * 3], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;
		//	TDOA_23 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER * 2], &data[CHANNEL_PER_BUFFER * 3], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;

		//	do_osh_algorithm(osh_azimuth, osh_elevation, TDOA_01, TDOA_20, TDOA_30, TDOA_12, TDOA_13, TDOA_23);
		//	cord3_trans(osh_azimuth, osh_elevation, s_theta_all, s_piangle_all, 20, mic_distance);

		//	do_AAI_algorithm2(aai_azimuth, aai_elevation, TDOA_01, TDOA_20, TDOA_30, TDOA_12, TDOA_13, TDOA_23);
		//	cord3_trans(aai_azimuth, aai_elevation, s_theta_all, s_piangle_all, 20, mic_distance);

		//	//fprintf(angle, "%lf\t%lf\t%lf\t%lf\n", osh_azimuth[4], osh_elevation[4], aai_azimuth[4], aai_elevation[4]);
		//	fprintf(stdout, "%lf\t%lf\t%lf\t%lf\n", osh_azimuth[4], osh_elevation[4], aai_azimuth[4], aai_elevation[4]);
		//	//fprintf(angle, "%lf\t%lf\t%lf\t%lf\t%lf\t", osh_azimuth[0], osh_azimuth[1], osh_azimuth[2], osh_azimuth[3], osh_azimuth[4]);
		//	//fprintf(angle, "%lf\t%lf\t%lf\t%lf\t%lf\t\t", aai_azimuth[0], aai_azimuth[1], aai_azimuth[2], aai_azimuth[3], aai_azimuth[4]);
		//	//fprintf(angle, "%lf\t%lf\t%lf\t%lf\t%lf\t", osh_elevation[0], osh_elevation[1], osh_elevation[2], osh_elevation[3], osh_elevation[4]);
		//	//fprintf(angle, "%lf\t%lf\t%lf\t%lf\t%lf\t\t", aai_elevation[0], aai_elevation[1], aai_elevation[2], aai_elevation[3], aai_elevation[4]);
		//	//fprintf(angle, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", TDOA_01, TDOA_12, TDOA_20, TDOA_30, TDOA_13, TDOA_23);

		//	i = 0;
		//	case_cnt++;
		//	getchar();
		//}
	}
	start = clock() - start;

	printf("ave : %lf %lf\n", error_azimuth / (double)case_cnt, error_elevation / (double)case_cnt);
	printf("\n%f seconds \n",((float)start) / CLOCKS_PER_SEC);

	fclose(record);
	fclose(angle);
	printf("\a\a\a");
}

#else READ_FROM_DAQ

#include <NIDAQmx.h>

// NI-DAQ 장비 세팅
#define DAQmxErrChk(functionCall) if( DAQmxFailed(error=(functionCall)) ) goto Error; else

// NI 콜백
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
		
	record = fopen("./data/2017-10-12/record/녹음_3m_a0도_e60도_25600hz.txt", "w");
	angle = fopen("./data/2017-10-12/angle/각도_3m_a0도_e60도_25600hz.txt", "w");

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
	float64     data[CHANNEL_COUNT * CHANNEL_PER_BUFFER] = { 0.0, };  // 4채널 데이터 버퍼

	// 샘플 딜레이, 타임 딜레이
	double TDOA_01 = 0, TDOA_20 = 0, TDOA_30 = 0, TDOA_12 = 0, TDOA_13 = 0, TDOA_23 = 0;

	// OSH, AAI 방위 고도
	double osh_azimuth[5] = { 0, }, osh_elevation[5] = { 0, };
	double aai_azimuth[5] = { 0, }, aai_elevation[5] = { 0, };
	double s_theta_all[10], s_piangle_all[10];

	// cpsp 인자
	double MAX_Delay = 40, temp = 999;

	// 알고리즘 성능 측정을 위한 변수
	static double case_cnt = 0, osh_error, aai_error = 0;
	static clock_t osh_start = 0, aai_start = 0, start = 0;
	static double osh_runtime = 0.0, aai_runtime = 0.0, runtime= 0.0;

	// DAQmx Read Code
	DAQmxErrChk(DAQmxReadAnalogF64(taskHandle, CHANNEL_PER_BUFFER, 10.0, DAQmx_Val_GroupByChannel, 
		data, CHANNEL_COUNT * CHANNEL_PER_BUFFER, &read, NULL));
	
	if (read >= CHANNEL_PER_BUFFER)
	{
		// printf("Acquired %d samples. Total %d\n", (int)read, (int)(totalRead += read));

		// 파일 저장
		for (i = 0; i < CHANNEL_PER_BUFFER; i++)
			fprintf(record, "%lf\t%lf\t%lf\t%lf\n", data[i], data[i + CHANNEL_PER_BUFFER], data[i + 2 * CHANNEL_PER_BUFFER], data[i + 3 * CHANNEL_PER_BUFFER]);

		start = clock();

		// 0. CPSP로 타임 딜레이 계산
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
		//osh_runtime += (double)(clock() - osh_start) / CLOCKS_PER_SEC; // 수행 시간 측정
	 
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

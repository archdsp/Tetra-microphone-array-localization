/*
작성자 : 최지수
작성일 : 2017-09-14
설명 : 정사면체 배치의 마이크로폰 어레이에 들어온 음원의 TDOA를 이용하여 음원의 방위각 및 고도각을 알아내는 코드로
AAI알고리즘의 일부와 오상헌 박사의 알고리즘을 혼용 및 배합하여 적용하였다.
아래 코드는 라인 당 6개의 TDOA가 있는 파일로부터 TDOA를 얻어 사용한다.
*/

#pragma once
#include <windows.h>
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

int main()
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

	double deg2rad = (M_PI / 180.0);
	double rad2deg = (180.0 / M_PI);

	double GT_azimuth = 0.0, GT_elevation = 0.0;
	int dist = 3, func_idx = 0;
	void(*algorithm[2])(double *planes_azimuth, double *planes_elevation,
		double TDOA_01, double TDOA_20, double TDOA_30, double TDOA_12, double TDOA_13, double TDOA_23);
	= { do_osh_algorithm, do_AAI_algorithm2 };
	
	char func_name[2][] = { "OSH", "AAI" };

	for (func_idx = 1; func_idx < 2; func_idx++)
	{
		for (dist = 3; dist < 9; dist += 3)
		{
			for (GT_elevation = 0.0; GT_elevation < 90; GT_elevation += 30.0)
			{
				for (GT_azimuth = 0.0; GT_azimuth < 360; GT_azimuth += 60.0)
				{
					double case_cnt = 1, detected = 0;
					int i = 0;

					// 주파수(1초당 증가되는 카운트수)를 구한다.
					/*LARGE_INTEGER liCounter1, liCounter2, liFrequency;
					QueryPerformanceFrequency(&liFrequency);
					double time = 0.0, time_sum = 0.0, time_avg = 0.0;*/

					double error_azimuth = 0.0, error_elevation = 0.0, err_tmp_azi = 0.0, err_tmp_ele;
					char file_name[] = "C:\\Users\\pblwo\\Desktop\\2017-10-14\\record\\고도각 00도\\녹음_3m_a000도_e00도_25600hz.txt";
					sprintf(file_name, "C:\\Users\\pblwo\\Desktop\\2017-10-14\\record\\고도각 %d도\\녹음_%dm_a%03d도_e%02d도_25600hz.txt", (int)GT_elevation, dist, (int)GT_azimuth, (int)GT_elevation);

					printf("file : %s\n", file_name);
					FILE* record = fopen(file_name, "r");
					if (record == NULL)
					{
						fprintf(stdout, "no file\n");
						return 0;
					}

					// 코드 수행...
					while (fscanf(record, "%lf\t%lf\t%lf\t%lf\n",
						&data[i], &data[i + CHANNEL_PER_BUFFER], &data[i + 2 * CHANNEL_PER_BUFFER], &data[i + 3 * CHANNEL_PER_BUFFER]) != EOF)
						//while(case_cnt <= 1000)
					{
						if (i == CHANNEL_PER_BUFFER)
						{
							//rand_TDOA(TDOA, mic_distance, &GT_azimuth, &GT_elevation);
							//TDOA_01 = TDOA[0];
							//TDOA_12 = TDOA[1];
							//TDOA_20 = TDOA[2];
							//TDOA_30 = TDOA[3];
							//TDOA_13 = TDOA[4];
							//TDOA_23 = TDOA[5];

							// 0. CPSP로 타임 딜레이 계산
							TDOA_01 = CPSP_FILT_normal(&data[0], &data[CHANNEL_PER_BUFFER], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;
							TDOA_12 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER], &data[CHANNEL_PER_BUFFER * 2], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;
							TDOA_20 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER * 2], &data[0], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;
							TDOA_30 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER * 3], &data[0], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;
							TDOA_13 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER], &data[CHANNEL_PER_BUFFER * 3], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;
							TDOA_23 = CPSP_FILT_normal(&data[CHANNEL_PER_BUFFER * 2], &data[CHANNEL_PER_BUFFER * 3], FFT_BUFF, MAX_Delay) / SAMPLING_RATE;

							// 코드 수행 전 카운트 저장
							//QueryPerformanceCounter(&liCounter1);

							algorithm[func_idx](azimuth, elevation, TDOA_01, TDOA_20, TDOA_30, TDOA_12, TDOA_13, TDOA_23);
							//do_osh_algorithm(azimuth, elevation, TDOA_01, TDOA_20, TDOA_30, TDOA_12, TDOA_13, TDOA_23);
							//do_AAI_algorithm2(azimuth, elevation, TDOA_01, TDOA_20, TDOA_30, TDOA_12, TDOA_13, TDOA_23);
							cord3_trans(azimuth, elevation, s_theta_all, s_piangle_all, 20, mic_distance);

							// 코드 수행 후 카운트 저장
							//QueryPerformanceCounter(&liCounter2);

							//time = (double)(liCounter2.QuadPart - liCounter1.QuadPart) / (double)liFrequency.QuadPart;
							//time_sum += time;

							// 수행 시간 계산값 출력
							//printf("평균 수행 시간 = %f 초\n", time_sum / case_cnt );

							//printf("%.0lf : %lf %lf\n", case_cnt, azimuth[4], elevation[4]);
							if (!isnan(azimuth[4]) || !isnan(elevation[4]))
							{
								azimuth[4] *= deg2rad;
								elevation[4] *= deg2rad;
								err_tmp_azi = fabs(((GT_azimuth * deg2rad) - azimuth[4]) * rad2deg);
								err_tmp_ele = fabs(((GT_elevation * deg2rad) - elevation[4]) * rad2deg);

								if (err_tmp_azi < 10 && err_tmp_ele < 10)
								{
									detected++;
									error_azimuth += err_tmp_azi;
									error_elevation += err_tmp_ele;
								}
							}

							case_cnt++;
							i = 0;
						}
						i++;
					} // frame while loop

					printf("detect err angle : %lf %lf\n", error_azimuth / detected, error_elevation / detected);
					printf("error rate : %.0lf / %.0lf, %lf\n\n", detected, case_cnt, detected / case_cnt * 100.0);

					if (record != NULL)
						fclose(record);
				}// var azimuth for loop
			}// var elevation for loop
		} // var dist for loop
	}
	
	printf("done\n");
	return 0;
}

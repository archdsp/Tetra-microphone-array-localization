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
#include <string.h>

#include "sort.h"
#include "CPSP.h"
#include "OSH_sound_localization.h"
#include "tetra_sound_localization.h"

// 채널 갯수
#define CHANNEL_COUNT 4

// 채널 당 데이터 갯수
#define CHANNEL_BUF 30000

// 음속
const double SOUND_SPEED = 340.0;

// 알고리즘 수행 시간
clock_t start = 0;
double runtime = 0.0;

int main(int argc, char* argv[])
{
	// 샘플링 주파수
	double samp_freq = 51200;

	// 마이크간 거리
	double mic_distance = 0.35;

	// 샘플 딜레이
	double Delay_01 = 0, Delay_20 = 0, Delay_12 = 0, Delay_30 = 0, Delay_13 = 0, Delay_31 = 0, Delay_23 = 0;
	
	// 타임 딜레이
	double TDOA_01 = 0, TDOA_20 = 0, TDOA_30 = 0, TDOA_12 = 0, TDOA_13 = 0, TDOA_23 = 0;

	// OSH 방위 고도
	double side_plane_based_azimuth[5] = { 0, }, side_plane_based_elevation[5] = { 0, };
	double s_theta_all[10], s_piangle_all[10];

	// AAI 방위 고도
	SoundLocation *location[CHANNEL_COUNT] = { NULL, };

	// AAI 알고리즘에서 쓰이는 첫번째와 두번째 마이크간 타임딜레이, 첫번째와 세번째 마이크간 타임딜레이 
	double t_F, t_L = 0.0;
	
	// 마이크 순서
	int mic_order = 0;

	// 알고리즘 성능 측정을 위한 변수
	double cnt = 1, detect_rate[2] = { 0.0, 0.0 }, error[4] = { 0.0, 0.0, 0.0, 0.0 };

	// 입력 음향 데이터 파일
	FILE *data_file = NULL;
	data_file = fopen("./data/integrated_data.txt", "r");

	start = clock();
	while (fscanf(data_file, "%lf %lf %lf %lf %lf %lf", &Delay_01, &Delay_12, &Delay_20, &Delay_30, &Delay_13, &Delay_23) != EOF)
	{
		// 타임 딜레이
		TDOA_01 = Delay_01 / samp_freq;
		TDOA_12 = Delay_12 / samp_freq;
		TDOA_20 = Delay_20 / samp_freq;
		TDOA_30 = Delay_30 / samp_freq;
		TDOA_13 = Delay_13 / samp_freq;
		TDOA_23 = Delay_23 / samp_freq;

		/* 각 평면에서 각도 구함 */

		/* OSH */
		//바닥면
		mic3_dir(&side_plane_based_azimuth[0], &side_plane_based_elevation[0], TDOA_01, TDOA_12, TDOA_20, mic_distance, SOUND_SPEED);		

		//옆면 A
		mic3_dir(&side_plane_based_azimuth[1], &side_plane_based_elevation[1], TDOA_30, TDOA_01, TDOA_13, mic_distance, SOUND_SPEED);

		//옆면 B
		mic3_dir(&side_plane_based_azimuth[2], &side_plane_based_elevation[2], -TDOA_13, TDOA_12, TDOA_23, mic_distance, SOUND_SPEED);

		//옆면 C
		mic3_dir(&side_plane_based_azimuth[3], &side_plane_based_elevation[3], -TDOA_23, TDOA_20, -TDOA_30, mic_distance, SOUND_SPEED);

		/* AAI */
		//// 바닥면
		//mic_order = find_tF_tL_from_cpsp(TDOA_30, TDOA_01, TDOA_13, &t_F, &t_L);
		//location[1] = AAI_direction(mic_distance, SOUND_SPEED, mic_order, t_F, t_L);
		//AAI2OSH(location[1], mic_distance);
		//
		//// 옆면 A
		//mic_order = find_tF_tL_from_cpsp(TDOA_01, TDOA_12, TDOA_20, &t_F, &t_L);
		//location[0] = AAI_direction(mic_distance, SOUND_SPEED, mic_order, t_F, t_L);
		//AAI2OSH(location[0], mic_distance);
		//
		//// 옆면 B
		//mic_order = find_tF_tL_from_cpsp(-TDOA_13, TDOA_12, TDOA_23, &t_F, &t_L);
		//location[2] = AAI_direction(mic_distance, SOUND_SPEED, mic_order, t_F, t_L);
		//AAI2OSH(location[2], mic_distance);
		//
		//// 옆면 C
		//mic_order = find_tF_tL_from_cpsp(-TDOA_23, TDOA_20, -TDOA_30, &t_F, &t_L);
		//location[3] = AAI_direction(mic_distance, SOUND_SPEED, mic_order, t_F, t_L);
		//AAI2OSH(location[3], mic_distance);
		//
		//// AAI로 구한 좌표를 대입 변환
		//side_plane_based_azimuth[0] = location[0]->azimuth;
		//side_plane_based_elevation[0] = location[0]->elevation;

		//side_plane_based_azimuth[1] = location[1]->azimuth;
		//side_plane_based_elevation[1] = location[1]->elevation;

		//side_plane_based_azimuth[2] = location[2]->azimuth;
		//side_plane_based_elevation[2] = location[2]->elevation;

		//side_plane_based_azimuth[3] = location[3]->azimuth;
		//side_plane_based_elevation[3] = location[3]->elevation;

		// 모든 면 바닥면 기준으로 좌표 변환
		cord3_trans(side_plane_based_azimuth, side_plane_based_elevation, s_theta_all, s_piangle_all, 20, mic_distance);

		printf("event = %lf degree\n", (cnt-1) * 30);	
		printf("azi = %lf ele = %lf\n", side_plane_based_azimuth[4], side_plane_based_elevation[4]);
		printf("=====================================\n\n");
		
		cnt++;
	}

	// 수행 시간 측정
	runtime = (double)(clock() - start) / CLOCKS_PER_SEC;
	printf("runtime : %lf\n", runtime);
	printf("average runtime : %lf\n", runtime / cnt);

	fclose(data_file);

	return 0;
}
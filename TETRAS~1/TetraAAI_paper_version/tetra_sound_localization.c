#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include  "tetra_sound_localization.h"

double SOUND_VEL = 340.0;
double mic_distance = 0.35;

void rand_TDOA(double TDOA[], double mic_dist, double *sound_azimuth, double* sound_elevation)
{
	double sound_dist = 0.0;
	double sound_x = 0.0, sound_y = 0.0, sound_z = 0.0;
	double tmp = 0.0;
	double TOA0 = 0.0, TOA1 = 0.0, TOA2 = 0.0, TOA3 = 0.0;
	double r = mic_dist / sqrt(3);
	double h = sqrt(6) / 3.0 * mic_distance;

	sound_dist = rand() % 10 + 1;
	*sound_azimuth = rand() % 360;
	*sound_elevation = rand() % 90;

	*sound_azimuth *= (M_PI / 180.0);
	*sound_elevation *= (M_PI / 180.0);

	sound_x = cos(*sound_azimuth);
	sound_y = sin(*sound_azimuth);
	tmp = sqrt(sound_x*sound_x + sound_y * sound_y);
	sound_z = tmp * tan(*sound_elevation);

	tmp = sqrt(sound_x*sound_x + sound_y * sound_y + sound_z*sound_z);
	sound_x = (sound_x / tmp) * sound_dist;
	sound_y = (sound_y / tmp) * sound_dist;
	sound_z = (sound_z / tmp) * sound_dist;

	TOA0 = sqrt(pow( sound_x - r, 2) + pow(sound_y - 0.0, 2) + pow(sound_z - 0.0, 2)) / SOUND_VEL;
	TOA1 = sqrt(pow(sound_x - (-r / 2.0), 2) + pow(sound_y - (mic_distance / 2.0), 2) + pow(sound_z - 0.0, 2)) / SOUND_VEL;
	TOA2 = sqrt(pow(sound_x - (-r / 2.0), 2) + pow(sound_y - (-mic_distance / 2.0), 2) + pow(sound_z - 0.0, 2)) / SOUND_VEL;
	TOA3 = sqrt(pow(sound_x - 0.0, 2) + pow(sound_y - 0.0, 2) + pow(sound_z - h, 2)) / SOUND_VEL;

	TDOA[0] = TOA1 - TOA0;
	TDOA[1] = TOA2 - TOA1;
	TDOA[2] = TOA0 - TOA2;
	TDOA[3] = TOA0 - TOA3;
	TDOA[4] = TOA3 - TOA1;
	TDOA[5] = TOA3 - TOA2;
}

void do_AAI_algorithm2(double *planes_azimuth, double *planes_elevation,
	double TDOA_01, double TDOA_20, double TDOA_30, double TDOA_12, double TDOA_13, double TDOA_23)
{
	// AAI 알고리즘에서 쓰이는 첫번째와 두번째 마이크간 타임딜레이, 첫번째와 세번째 마이크간 타임딜레이 
	double t_F, t_L = 0.0;

	// 마이크  입력 순서
	int mic_order = 0;

	mic_order = find_tF_tL_from_cpsp(TDOA_01, TDOA_12, TDOA_20, &t_F, &t_L);
	AAI_direction2(&planes_azimuth[0], &planes_elevation[0], mic_distance, SOUND_VEL, mic_order, t_F, t_L);

	mic_order = find_tF_tL_from_cpsp(TDOA_30, TDOA_01, TDOA_13, &t_F, &t_L);
	AAI_direction2(&planes_azimuth[1], &planes_elevation[1], mic_distance, SOUND_VEL, mic_order, t_F, t_L);

	mic_order = find_tF_tL_from_cpsp(-TDOA_13, TDOA_12, TDOA_23, &t_F, &t_L);
	AAI_direction2(&planes_azimuth[2], &planes_elevation[2], mic_distance, SOUND_VEL, mic_order, t_F, t_L);

	mic_order = find_tF_tL_from_cpsp(-TDOA_23, TDOA_20, -TDOA_30, &t_F, &t_L);
	AAI_direction2(&planes_azimuth[3], &planes_elevation[3], mic_distance, SOUND_VEL, mic_order, t_F, t_L);
}

void do_osh_algorithm(double *planes_azimuth, double *planes_elevation, 
	double TDOA_01, double TDOA_20, double TDOA_30, double TDOA_12, double TDOA_13, double TDOA_23)
{
	/* OSH : 바닥면, 옆면 A, B, C */
	mic3_dir(&planes_azimuth[0], &planes_elevation[0], TDOA_01, TDOA_12, TDOA_20, mic_distance, SOUND_VEL);
	mic3_dir(&planes_azimuth[1], &planes_elevation[1], TDOA_30, TDOA_01, TDOA_13, mic_distance, SOUND_VEL);
	mic3_dir(&planes_azimuth[2], &planes_elevation[2], -TDOA_13, TDOA_12, TDOA_23, mic_distance, SOUND_VEL);
	mic3_dir(&planes_azimuth[3], &planes_elevation[3], -TDOA_23, TDOA_20, -TDOA_30, mic_distance, SOUND_VEL);
}

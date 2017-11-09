#pragma once

#include "OSH_sound_localization.h"
#include "AAI_sound_localization.h"

#define CHANNEL_COUNT 4


void rand_TDOA(double TDOA[], double mic_dist, double *sound_azimuth, double* sound_elevation);

// 음속
extern double SOUND_VEL;

// 마이크간 거리
extern double mic_distance;

void do_AAI_algorithm2(double *planes_azimuth, double *planes_elevation,
	double TDOA_01, double TDOA_20, double TDOA_30, double TDOA_12, double TDOA_13, double TDOA_23);

// OSH 알고리즘
void do_osh_algorithm(double *planes_azimuth, double *planes_elevation,
	double TDOA_01, double TDOA_20, double TDOA_30, double TDOA_12, double TDOA_13, double TDOA_23);

// AAI 알고리즘
void do_AAI_algorithm(double *planes_azimuth, double *planes_elevation,
	double TDOA_01, double TDOA_20, double TDOA_30, double TDOA_12, double TDOA_13, double TDOA_23);
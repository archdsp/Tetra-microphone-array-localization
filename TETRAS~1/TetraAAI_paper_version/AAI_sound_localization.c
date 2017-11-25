#pragma once

#include "AAI_sound_localization.h"

void AAI_direction2(double *planes_azimuth, double *planes_elevation, double S, double sound_speed, int mic_order, double t_F, double t_L)
{
	double S2 = t_F * sound_speed;
	double S3 = t_L * sound_speed;
	double x = 0.0, y = 0.0, z = 0.0;
	int i = 0;
	double rad = 0, degree = 0.0;
	double tmp_location[3] = { 0.0, 0.0, 0.0 };

	x = S3 * (S2 - S3) / S;
	y = -(S2*S3 + pow(S3, 2)) / (S * sqrt(3));
	double z_tmp = pow(S3, 2)*(pow(S, 2) - pow(S3 - S2, 2)) / pow(S, 2) - pow(y, 2);
	z = sqrt(fabs(z_tmp));

	// (0 ~ 60 degree)
	if (mic_order == MIC_ORDER_012)
	{
		x *= -1;
		degree = 90.0;
	}

	// (0 ~ 120 degree)
	else if (mic_order == MIC_ORDER_102)
		degree = -150.0;

	// (120 ~ 180 degree)
	else if (mic_order == MIC_ORDER_120)
	{
		x *= -1;
		degree = -150.0;
	}

	// (180 ~ 240 degree)
	else if (mic_order == MIC_ORDER_210)
		degree = -30.0;

	// (240 ~ 300 degree)
	else if (mic_order == MIC_ORDER_201)
	{
		x *= -1;
		degree = -30.0;
	}

	// (300 ~ 360 degree)
	else if (mic_order == MIC_ORDER_021)
		degree = 90.0;

	rad = degree * M_PI / 180.0;

	// 회전변환 행렬
	double rot_mat[2][2] = { { cos(rad), -sin(rad)},
	{ sin(rad), cos(rad)} };

	// 회전변환 수행
	for (i = 0; i < 2; i++)
		tmp_location[i] = rot_mat[i][0] * x + rot_mat[i][1] * y;

	x = tmp_location[0];
	y = tmp_location[1];
//	z = tmp_location[2];

	*planes_azimuth = atan2(y, x) * 180.0 / M_PI;
	if (*planes_azimuth < 0)
		*planes_azimuth += 360.0;

	if (*planes_azimuth == 360)
		*planes_azimuth = 0;

	*planes_elevation = atan2(z, sqrt(pow(x, 2) + pow(y, 2))) * 180.0 / M_PI;
}


// 마이크 도착 순서, t_F, t_L 찾음
// @param t_F : 첫번째로 소리를 받은 센서와 두번째로 소리를 받은 센서간  시간 차
// @param t_L : 첫번째로 소리를 받은 센서와 마지막로 소리를 받은 센서간  시간 차
// @return : 마이크 도착 순서
int find_tF_tL_from_cpsp(double cpsp_base_TDOA01, double cpsp_base_TDOA12, double cpsp_base_TDOA20, double *t_F, double * t_L)
{
	// mic arrival sequence : 0->1->2 (0 ~ 60 degree)
	if (cpsp_base_TDOA01 >= 0 && cpsp_base_TDOA12 >= 0) 
	{
		*t_F = fabs(cpsp_base_TDOA01);
		*t_L = fabs(cpsp_base_TDOA20); 
		return MIC_ORDER_012;
	}
	
	// mic arrival sequence : 1->0->2 (60 ~ 120 degree)
	else if (cpsp_base_TDOA01 <= 0 && cpsp_base_TDOA20 <= 0)
	{
		*t_F = fabs(cpsp_base_TDOA01);
		*t_L = fabs(cpsp_base_TDOA12);
		return MIC_ORDER_102;
	}
	
	// mic arrival sequence : 1->2->0 (120 ~ 180 degree)
	else if (cpsp_base_TDOA12 >= 0 && cpsp_base_TDOA20 >= 0)
	{
		*t_F = fabs(cpsp_base_TDOA12);
		*t_L = fabs(cpsp_base_TDOA01);
		return MIC_ORDER_120;
	}

	// mic arrival sequence : 2->1->0 (180 ~ 240 degree)
	else if (cpsp_base_TDOA01 <= 0 && cpsp_base_TDOA12 <= 0)
	{
		*t_F = fabs(cpsp_base_TDOA12);
		*t_L = fabs(cpsp_base_TDOA20);
		return MIC_ORDER_210;
	}
	
	// mic arrival sequence : 2->0->1 (240 ~ 300 degree)
	else if (cpsp_base_TDOA01 >= 0 && cpsp_base_TDOA20 >= 0)
	{
		*t_F = fabs(cpsp_base_TDOA20);
		*t_L = fabs(cpsp_base_TDOA12);
		return MIC_ORDER_201;
	}

	// mic arrival sequence : 0->2->1 (300 ~ 360 degree)
	else if (cpsp_base_TDOA12 <= 0 && cpsp_base_TDOA20 <= 0)
	{
		*t_F = fabs(cpsp_base_TDOA20);
		*t_L = fabs(cpsp_base_TDOA01);
		return MIC_ORDER_021;
	}
}

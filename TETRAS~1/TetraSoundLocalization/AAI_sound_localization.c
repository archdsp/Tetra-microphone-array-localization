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

// AAI 알고리즘으로 방향각, 고도각, 좌표 알아내기
// @param S : 센서사이 거리 
// @param samp_freq : 샘플링 레이트
// @param sound_speed : 음속
// @param t_F : 첫번째로 소리를 받은 센서와 두번째로 소리를 받은 센서간  시간 차
// @param t_L : 첫번째로 소리를 받은 센서와 마지막로 소리를 받은 센서간  시간 차
// 메모리 할당하여 SoundLocation* 로 리턴
SoundLocation* AAI_direction(double S, double sound_speed, int mic_order, double t_F, double t_L)
{
	SoundLocation* ret = NULL;
	ret = (SoundLocation*)malloc(sizeof(SoundLocation));

	double S2 = t_F * sound_speed;
	double S3 = t_L * sound_speed;

	ret->mic_order = mic_order;

	ret->x = S3 * (S2 - S3) / S;
	ret->y = -(S2*S3 + pow(S3, 2)) / (S * sqrt(3));
	double z_tmp = pow(S3, 2)*(pow(S, 2) - pow(S3 - S2, 2)) / pow(S, 2) - pow(ret->y, 2);
	ret->z = sqrt(fabs(z_tmp));
	
	ret->azimuth = atan2(ret->y, ret->x) * 180.0 / M_PI;
	
	if (ret->azimuth < 0) 
		ret->azimuth += 360.0;
	
	ret->elevation = atan(ret->z / sqrt(pow(ret->x, 2) + pow(ret->y, 2))) * 180 / M_PI;

	return ret;
}

SoundLocation* AAI2OSH(SoundLocation* AAI_based_location)
{	
	// (0 ~ 60 degree)
	if (AAI_based_location->mic_order == MIC_ORDER_012)
	{
		AAI_based_location->x *= -1;
		location_3d_rotation(AAI_based_location, 90.0);
	}

	// (0 ~ 120 degree)
	else if (AAI_based_location->mic_order == MIC_ORDER_102)
		location_3d_rotation(AAI_based_location, -150.0);

	// (120 ~ 180 degree)
	else if (AAI_based_location->mic_order == MIC_ORDER_120)
	{
		AAI_based_location->x *= -1;
		location_3d_rotation(AAI_based_location, -150.0);
	}

	// (180 ~ 240 degree)
	else if (AAI_based_location->mic_order == MIC_ORDER_210)
		location_3d_rotation(AAI_based_location, -30.0);

	// (240 ~ 300 degree)
	else if (AAI_based_location->mic_order == MIC_ORDER_201)
	{
		AAI_based_location->x *= -1;
		location_3d_rotation(AAI_based_location, -30.0);
	}

	// (300 ~ 360 degree)
	else if (AAI_based_location->mic_order == MIC_ORDER_021)
		location_3d_rotation(AAI_based_location, 90.0);

	else
		return;

	calc_azimuth_elevation(AAI_based_location);
	
	return AAI_based_location;
}

// 방향각, 고도각 계산
void calc_azimuth_elevation(SoundLocation* location)
{
	location->azimuth = atan2(location->y, location->x) * 180.0 / M_PI;
	if (location->azimuth < 0)
		location->azimuth += 360.0;

	if (location->azimuth == 360)
		location->azimuth = 0;

	location->elevation = atan2(location->z, sqrt(pow(location->x, 2) + pow(location->y, 2))) * 180.0 / M_PI;
}

// 3 차원 회전 변환
void location_3d_rotation(SoundLocation* location, double degree)
{
	int i = 0;
	double rad = degree*M_PI / 180.0;
	double tmp_location[3] = {0.};

	// 회전변환 행렬
	double rot_mat[3][3] = {{cos(rad), -sin(rad),  0},
							{sin(rad), cos(rad),   0}, 
							{0,        0,          1} };

	// 회전변환 수행
	for (i = 0; i <  3; i ++)
		tmp_location[i] = rot_mat[i][0] * location->x + rot_mat[i][1] * location->y + rot_mat[i][2] * location->z;
	
	location->x = tmp_location[0];
	location->y = tmp_location[1];
	location->z = tmp_location[2];
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

#pragma once
#include <malloc.h>
#include <math.h>

#define M_PI       3.14159265358979323846

#define MIC_ORDER_012 12
#define MIC_ORDER_021 21
#define MIC_ORDER_102 102
#define MIC_ORDER_120 120
#define MIC_ORDER_201 201
#define MIC_ORDER_210 210


// AAI 알고리즘으로 방향각, 고도각, 좌표 알아내기
// @param S : 센서사이 거리 
// @param sound_speed : 음속
// @param mic_order : 마이크 순서
// @param t_F : 첫번째로 소리를 받은 센서와 두번째로 소리를 받은 센서간  시간 차
// @param t_L : 첫번째로 소리를 받은 센서와 마지막로 소리를 받은 센서간  시간 차
void AAI_direction2(double *planes_azimuth, double *planes_elevation, double S, double sound_speed, int mic_order, double t_F, double t_L);

// 마이크 도착 순서, t_F, t_L 찾음
// @param t_F : 첫번째로 소리를 받은 센서와 두번째로 소리를 받은 센서간  시간 차
// @param t_L : 첫번째로 소리를 받은 센서와 마지막로 소리를 받은 센서간  시간 차
// @return : 마이크 도착 순서
int find_tF_tL_from_cpsp(double cpsp_base_TDOA01, double cpsp_base_TDOA12, double cpsp_base_TDOA20, double *t_F, double * t_L);
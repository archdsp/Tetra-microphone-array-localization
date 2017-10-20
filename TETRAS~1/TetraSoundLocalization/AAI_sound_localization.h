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

// ���� �ҽ� ��ġ
typedef struct _SoundLocation
{
	int mic_order;
	double azimuth;
	double elevation;
	double x;
	double y;
	double z;
}SoundLocation;


void AAI_direction2(double *planes_azimuth, double *planes_elevation, double S, double sound_speed, int mic_order, double t_F, double t_L);

// AAI �˰������� ���Ⱒ, ����, ��ǥ �˾Ƴ���
// @param S : �������� �Ÿ� 
// @param sound_speed : ����
// @param mic_order : ����ũ ����
// @param t_F : ù��°�� �Ҹ��� ���� ������ �ι�°�� �Ҹ��� ���� ������  �ð� ��
// @param t_L : ù��°�� �Ҹ��� ���� ������ �������� �Ҹ��� ���� ������  �ð� ��
SoundLocation* AAI_direction(double S, double sound_speed, int mic_order, double t_F, double t_L);

// AAI ��� ��ǥ���� ������ �ڻ� �˰������� ����
SoundLocation* AAI2OSH(SoundLocation* AAI_based_location);

// ���Ⱒ, ���� ���
void calc_azimuth_elevation(SoundLocation* location);

// ȸ����ȯ
void location_3d_rotation(SoundLocation* location, double degree);

// ����ũ ���� ����, t_F, t_L ã��
// @param t_F : ù��°�� �Ҹ��� ���� ������ �ι�°�� �Ҹ��� ���� ������  �ð� ��
// @param t_L : ù��°�� �Ҹ��� ���� ������ �������� �Ҹ��� ���� ������  �ð� ��
// @return : ����ũ ���� ����
int find_tF_tL_from_cpsp(double cpsp_base_TDOA01, double cpsp_base_TDOA12, double cpsp_base_TDOA20, double *t_F, double * t_L);
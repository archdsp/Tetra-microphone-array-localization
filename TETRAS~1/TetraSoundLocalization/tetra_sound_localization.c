#include "tetra_sound_localization.h"
#include <stdio.h>
double SOUND_VEL = 340.0;
double mic_distance = 0.35;

void do_AAI_algorithm2(double *planes_azimuth, double *planes_elevation,
	double TDOA_01, double TDOA_20, double TDOA_30, double TDOA_12, double TDOA_13, double TDOA_23)
{
	// AAI �˰��򿡼� ���̴� ù��°�� �ι�° ����ũ�� Ÿ�ӵ�����, ù��°�� ����° ����ũ�� Ÿ�ӵ����� 
	double t_F, t_L = 0.0;

	// ����ũ  �Է� ����
	int mic_order = 0;

	mic_order = find_tF_tL_from_cpsp(TDOA_01, TDOA_12, TDOA_20, &t_F, &t_L);
	AAI_direction2(&planes_azimuth[0], &planes_elevation[0], mic_distance, SOUND_VEL, mic_order, t_F, t_L);

	mic_order = find_tF_tL_from_cpsp(TDOA_30, TDOA_01, TDOA_13, &t_F, &t_L);
	AAI_direction2(&planes_azimuth[1], &planes_elevation[1], mic_distance, SOUND_VEL, mic_order, t_F, t_L);

	mic_order = find_tF_tL_from_cpsp(-TDOA_13, TDOA_12, TDOA_23, &t_F, &t_L);
	AAI_direction2(&planes_azimuth[2], &planes_elevation[2], mic_distance, SOUND_VEL, mic_order, t_F, t_L);

	mic_order = find_tF_tL_from_cpsp(-TDOA_23, TDOA_20, -TDOA_30, &t_F, &t_L);
	AAI_direction2(&planes_azimuth[3], &planes_elevation[3], mic_distance, SOUND_VEL, mic_order, t_F, t_L);


	fprintf(stdout, "AAI : %lf\t%lf\t%lf\t%lf\n", planes_azimuth[0], planes_azimuth[1], planes_azimuth[2], planes_azimuth[3]);
}

void do_osh_algorithm(double *planes_azimuth, double *planes_elevation, 
	double TDOA_01, double TDOA_20, double TDOA_30, double TDOA_12, double TDOA_13, double TDOA_23)
{
	/* OSH : �ٴڸ�, ���� A, B, C */
	mic3_dir(&planes_azimuth[0], &planes_elevation[0], TDOA_01, TDOA_12, TDOA_20, mic_distance, SOUND_VEL);
	mic3_dir(&planes_azimuth[1], &planes_elevation[1], TDOA_30, TDOA_01, TDOA_13, mic_distance, SOUND_VEL);
	mic3_dir(&planes_azimuth[2], &planes_elevation[2], -TDOA_13, TDOA_12, TDOA_23, mic_distance, SOUND_VEL);
	mic3_dir(&planes_azimuth[3], &planes_elevation[3], -TDOA_23, TDOA_20, -TDOA_30, mic_distance, SOUND_VEL);
	fprintf(stdout, "OSH : %lf\t%lf\t%lf\t%lf\n", planes_azimuth[0], planes_azimuth[1], planes_azimuth[2], planes_azimuth[3]);
}

void do_AAI_algorithm(double *planes_azimuth, double *planes_elevation,
	double TDOA_01, double TDOA_20, double TDOA_30, double TDOA_12, double TDOA_13, double TDOA_23)
{
	// AAI �˰��򿡼� ���̴� ù��°�� �ι�° ����ũ�� Ÿ�ӵ�����, ù��°�� ����° ����ũ�� Ÿ�ӵ����� 
	double t_F, t_L = 0.0;

	// ����ũ  �Է� ����
	int mic_order = 0;

	// AAI ��� ���� �� ����
	SoundLocation *location[CHANNEL_COUNT] = { NULL, };

	/* AAI : �ٴڸ�, ���� A, B, C */
	mic_order = find_tF_tL_from_cpsp(TDOA_01, TDOA_12, TDOA_20, &t_F, &t_L);
	location[0] = AAI_direction(mic_distance, SOUND_VEL, mic_order, t_F, t_L);
	AAI2OSH(location[0]);

	mic_order = find_tF_tL_from_cpsp(TDOA_30, TDOA_01, TDOA_13, &t_F, &t_L);
	location[1] = AAI_direction(mic_distance, SOUND_VEL, mic_order, t_F, t_L);
	AAI2OSH(location[1]);

	mic_order = find_tF_tL_from_cpsp(-TDOA_13, TDOA_12, TDOA_23, &t_F, &t_L);
	location[2] = AAI_direction(mic_distance, SOUND_VEL, mic_order, t_F, t_L);
	AAI2OSH(location[2]);

	mic_order = find_tF_tL_from_cpsp(-TDOA_23, TDOA_20, -TDOA_30, &t_F, &t_L);
	location[3] = AAI_direction(mic_distance, SOUND_VEL, mic_order, t_F, t_L);
	AAI2OSH(location[3]);

	// AAI�� ���� ��ǥ�� ���� ��ȯ
	planes_azimuth[0] = location[0]->azimuth;
	planes_elevation[0] = location[0]->elevation;

	planes_azimuth[1] = location[1]->azimuth;
	planes_elevation[1] = location[1]->elevation;

	planes_azimuth[2] = location[2]->azimuth;
	planes_elevation[2] = location[2]->elevation;

	planes_azimuth[3] = location[3]->azimuth;
	planes_elevation[3] = location[3]->elevation;

}

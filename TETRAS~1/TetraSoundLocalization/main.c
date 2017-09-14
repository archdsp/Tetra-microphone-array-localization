/*
	�ۼ��� : ������
	�ۼ��� : 2017-09-14
	���� : �����ü ��ġ�� ����ũ���� ��̿� ���� ������ TDOA�� �̿��Ͽ� ������ ������ �� ������ �˾Ƴ��� �ڵ��
	AAI�˰����� �Ϻο� ������ �ڻ��� �˰����� ȥ�� �� �����Ͽ� �����Ͽ���.
	�Ʒ� �ڵ�� ���� �� 6���� TDOA�� �ִ� ���Ϸκ��� TDOA�� ��� ����Ѵ�.
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

// ä�� ����
#define CHANNEL_COUNT 4

// ä�� �� ������ ����
#define CHANNEL_BUF 30000

// ����
const double SOUND_SPEED = 340.0;

// �˰��� ���� �ð�
clock_t start = 0;
double runtime = 0.0;

int main(int argc, char* argv[])
{
	// ���ø� ���ļ�
	double samp_freq = 51200;

	// ����ũ�� �Ÿ�
	double mic_distance = 0.35;

	// ���� ������
	double Delay_01 = 0, Delay_20 = 0, Delay_12 = 0, Delay_30 = 0, Delay_13 = 0, Delay_31 = 0, Delay_23 = 0;
	
	// Ÿ�� ������
	double TDOA_01 = 0, TDOA_20 = 0, TDOA_30 = 0, TDOA_12 = 0, TDOA_13 = 0, TDOA_23 = 0;

	// OSH ���� ��
	double side_plane_based_azimuth[5] = { 0, }, side_plane_based_elevation[5] = { 0, };
	double s_theta_all[10], s_piangle_all[10];

	// AAI ���� ��
	SoundLocation *location[CHANNEL_COUNT] = { NULL, };

	// AAI �˰��򿡼� ���̴� ù��°�� �ι�° ����ũ�� Ÿ�ӵ�����, ù��°�� ����° ����ũ�� Ÿ�ӵ����� 
	double t_F, t_L = 0.0;
	
	// ����ũ ����
	int mic_order = 0;

	// �˰��� ���� ������ ���� ����
	double cnt = 1, detect_rate[2] = { 0.0, 0.0 }, error[4] = { 0.0, 0.0, 0.0, 0.0 };

	// �Է� ���� ������ ����
	FILE *data_file = NULL;
	data_file = fopen("./data/integrated_data.txt", "r");

	start = clock();
	while (fscanf(data_file, "%lf %lf %lf %lf %lf %lf", &Delay_01, &Delay_12, &Delay_20, &Delay_30, &Delay_13, &Delay_23) != EOF)
	{
		// Ÿ�� ������
		TDOA_01 = Delay_01 / samp_freq;
		TDOA_12 = Delay_12 / samp_freq;
		TDOA_20 = Delay_20 / samp_freq;
		TDOA_30 = Delay_30 / samp_freq;
		TDOA_13 = Delay_13 / samp_freq;
		TDOA_23 = Delay_23 / samp_freq;

		/* �� ��鿡�� ���� ���� */

		/* OSH */
		//�ٴڸ�
		mic3_dir(&side_plane_based_azimuth[0], &side_plane_based_elevation[0], TDOA_01, TDOA_12, TDOA_20, mic_distance, SOUND_SPEED);		

		//���� A
		mic3_dir(&side_plane_based_azimuth[1], &side_plane_based_elevation[1], TDOA_30, TDOA_01, TDOA_13, mic_distance, SOUND_SPEED);

		//���� B
		mic3_dir(&side_plane_based_azimuth[2], &side_plane_based_elevation[2], -TDOA_13, TDOA_12, TDOA_23, mic_distance, SOUND_SPEED);

		//���� C
		mic3_dir(&side_plane_based_azimuth[3], &side_plane_based_elevation[3], -TDOA_23, TDOA_20, -TDOA_30, mic_distance, SOUND_SPEED);

		/* AAI */
		//// �ٴڸ�
		//mic_order = find_tF_tL_from_cpsp(TDOA_30, TDOA_01, TDOA_13, &t_F, &t_L);
		//location[1] = AAI_direction(mic_distance, SOUND_SPEED, mic_order, t_F, t_L);
		//AAI2OSH(location[1], mic_distance);
		//
		//// ���� A
		//mic_order = find_tF_tL_from_cpsp(TDOA_01, TDOA_12, TDOA_20, &t_F, &t_L);
		//location[0] = AAI_direction(mic_distance, SOUND_SPEED, mic_order, t_F, t_L);
		//AAI2OSH(location[0], mic_distance);
		//
		//// ���� B
		//mic_order = find_tF_tL_from_cpsp(-TDOA_13, TDOA_12, TDOA_23, &t_F, &t_L);
		//location[2] = AAI_direction(mic_distance, SOUND_SPEED, mic_order, t_F, t_L);
		//AAI2OSH(location[2], mic_distance);
		//
		//// ���� C
		//mic_order = find_tF_tL_from_cpsp(-TDOA_23, TDOA_20, -TDOA_30, &t_F, &t_L);
		//location[3] = AAI_direction(mic_distance, SOUND_SPEED, mic_order, t_F, t_L);
		//AAI2OSH(location[3], mic_distance);
		//
		//// AAI�� ���� ��ǥ�� ���� ��ȯ
		//side_plane_based_azimuth[0] = location[0]->azimuth;
		//side_plane_based_elevation[0] = location[0]->elevation;

		//side_plane_based_azimuth[1] = location[1]->azimuth;
		//side_plane_based_elevation[1] = location[1]->elevation;

		//side_plane_based_azimuth[2] = location[2]->azimuth;
		//side_plane_based_elevation[2] = location[2]->elevation;

		//side_plane_based_azimuth[3] = location[3]->azimuth;
		//side_plane_based_elevation[3] = location[3]->elevation;

		// ��� �� �ٴڸ� �������� ��ǥ ��ȯ
		cord3_trans(side_plane_based_azimuth, side_plane_based_elevation, s_theta_all, s_piangle_all, 20, mic_distance);

		printf("event = %lf degree\n", (cnt-1) * 30);	
		printf("azi = %lf ele = %lf\n", side_plane_based_azimuth[4], side_plane_based_elevation[4]);
		printf("=====================================\n\n");
		
		cnt++;
	}

	// ���� �ð� ����
	runtime = (double)(clock() - start) / CLOCKS_PER_SEC;
	printf("runtime : %lf\n", runtime);
	printf("average runtime : %lf\n", runtime / cnt);

	fclose(data_file);

	return 0;
}
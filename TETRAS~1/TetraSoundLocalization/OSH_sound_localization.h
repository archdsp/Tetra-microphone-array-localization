#pragma once

#ifndef M_PI
	#define M_PI       3.14159265358979323846
#endif // !M_PI


// �� �迭���� ��ġ�� ������ ����� ����
double intsec_mean(double *A, double *B, int len);
double mic3_dir(double *theta, double *piangle, double Td_cl, double Td_lr, double Td_rc, double mic_distance, double sound_speed);
double cord3_trans(double *side_plane_based_azimuth, double *side_plane_based_elevation, double *s_theta_all, double *s_piangle_all, int angle_range, double mic_distance);
#pragma once
#include <math.h>
#include <stdio.h>
#include "OSH_sound_localization.h"
#include "sort.h"

// 두 배열에서 겹치는 값들의 평균을 구함
double intsec_mean(double *A, double *B, int len)
{
	int i, j, cnt;
	double C[50], avg = 0.0;

	cnt = 0;
	for (i = 0; i < len; i++)
		for (j = 0; j < len; j++)
			if (A[i] == B[j])
			{
				C[cnt] = A[i];
				cnt++;
			}

	if (cnt == 0)
		return NO_INTERSECT_VALUE;

	for (i = 0; i < cnt; i++)
		avg = avg + C[i];
	avg = avg / (double)cnt;

	return avg;
}

double mic3_dir(double *theta, double *piangle, double Td_cl, double Td_lr, double Td_rc, double mic_distance, double sound_speed)
{
	int loop = 0, loop2 = 0, loop3 = 0;
	double r, p_i, temp_arg;

	double a, A[4] = { 0, }, AA[6] = { 0, };
	double tempA[2] = { 0, }, tempA2[6] = { 0, };
	int Data_Acnt = 0, Data_Acnt2 = 0;
	int AAb[180] = { 0, };

	double b, B[4] = { 0, }, BB[6] = { 0,0,0,0,0,0 };
	double tempB[2] = { 0,0 }, tempB2[6] = { 0,0,0,0 };
	int Data_Bcnt = 0, Data_Bcnt2 = 0;
	int BBb[180] = { 0, };
	double c, C[4] = { 0, }, CC[6] = { 0, };
	double tempC[2] = { 0, }, tempC2[6] = { 0, };
	double temp, temp2;
	int loop_temp = 0, loop_temp2 = 0;
	int Data_Ccnt = 0, Data_Ccnt2 = 0;
	int CCb[180] = { 0, };

	int theta_temp1[180] = { 0, };
	int theta_diff1[180] = { 0, };

	int theta_diff1_idx[10] = { 0, };
	int theta_array[10] = { 0, };

	int theta_temp1_cnt = 0, AAb_cnt = 0, BBb_cnt = 0, CCb_cnt = 0;

	r = mic_distance / sqrt(3.0);

	//%9-1 식
	temp_arg = (sqrt(2.0)*sound_speed / (3 * r)) * sqrt(((Td_cl)*(Td_cl)) + ((Td_lr)*(Td_lr)) + ((Td_rc)*(Td_rc)));
	if (temp_arg>1) 
		temp_arg = 1;
	
	else if (temp_arg<-1) 
		temp_arg = -1;
	
	p_i = acos(temp_arg);

	//9-2식
	temp_arg = sound_speed*Td_cl / (sqrt(3.0)*r*cos(p_i));
	if (temp_arg>1) 
		temp_arg = 1;
	
	else if (temp_arg<-1) 
		temp_arg = -1;
	
	a = acos(temp_arg);

	A[0] = a - (M_PI / 6);
	A[1] = -a - (M_PI / 6);

	tempA[0] = (A[0] / M_PI) * 180;
	tempA[1] = (A[1] / M_PI) * 180;

	//AA=[AA AA+360 AA-360];
	tempA2[0] = tempA[0];
	tempA2[1] = tempA[1];
	tempA2[2] = tempA[0] + 360;
	tempA2[3] = tempA[1] + 360;
	tempA2[4] = tempA[0] - 360;
	tempA2[5] = tempA[1] - 360;

	//AA=AA(find(AA>=-360&AA<=360));
	Data_Acnt = 0;
	for (loop = 0; loop < 6; loop++)
		if (fabs(tempA2[loop]) <= 360)
		{
			AA[Data_Acnt] = tempA2[loop];
			Data_Acnt++;
		}

	//9-2식
	temp_arg = sound_speed*Td_lr / (sqrt(3.0)*r*cos(p_i));
	if (temp_arg>1)
		temp_arg = 1;

	else if (temp_arg<-1)
		temp_arg = -1;

	b = acos(temp_arg);

	B[0] = b + (M_PI / 2);
	B[1] = -b + (M_PI / 2);

	tempB[0] = (B[0] / M_PI) * 180;
	tempB[1] = (B[1] / M_PI) * 180;

	//AA=[AA AA+360 AA-360];
	tempB2[0] = tempB[0];
	tempB2[1] = tempB[1];
	tempB2[2] = tempB[0] + 360;
	tempB2[3] = tempB[1] + 360;
	tempB2[4] = tempB[0] - 360;
	tempB2[5] = tempB[1] - 360;

	//AA=AA(find(AA>=-360&AA<=360));
	Data_Bcnt = 0;
	for (loop = 0; loop < 6; loop++)
		if (fabs(tempB2[loop]) <= 360)
		{
			BB[Data_Bcnt] = tempB2[loop];
			Data_Bcnt++;
		}

	temp_arg = sound_speed*Td_rc / (sqrt(3.0)*r*cos(p_i));
	if (temp_arg>1)
		temp_arg = 1;

	else if (temp_arg<-1)
		temp_arg = -1;

	c = acos(temp_arg);

	C[0] = c - ((M_PI * 5) / 6);
	C[1] = -c - ((M_PI * 5) / 6);

	tempC[0] = (C[0] / M_PI) * 180;
	tempC[1] = (C[1] / M_PI) * 180;

	//AA=[AA AA+360 AA-360];
	tempC2[0] = tempC[0];
	tempC2[1] = tempC[1];
	tempC2[2] = tempC[0] + 360;
	tempC2[3] = tempC[1] + 360;
	tempC2[4] = tempC[0] - 360;
	tempC2[5] = tempC[1] - 360;

	//AA=AA(find(AA>=-360&AA<=360));
	Data_Ccnt = 0;
	for (loop = 0; loop < 6; loop++)
		if (fabs(tempC2[loop]) <= 360)
		{
			CC[Data_Ccnt] = tempC2[loop];
			Data_Ccnt++;
		}

	if (Data_Acnt > 0)
		Data_Acnt2 = Data_Acnt - 1;

	if (Data_Bcnt > 0)
		Data_Bcnt2 = Data_Bcnt - 1;

	if (Data_Ccnt > 0)
		Data_Ccnt2 = Data_Ccnt - 1;

	for (loop = -15; loop <= 15; loop++)
	{
		for (loop2 = 0; loop2 < Data_Acnt; loop2++)
			AAb[loop2 + ((loop + 15)*Data_Acnt)] = AA[loop2] + loop;

		for (loop2 = 0; loop2 < Data_Bcnt; loop2++)
			BBb[loop2 + ((loop + 15)*Data_Bcnt)] = BB[loop2] + loop;

		for (loop2 = 0; loop2 < Data_Ccnt; loop2++)
			CCb[loop2 + ((loop + 15)*Data_Ccnt)] = CC[loop2] + loop;

	}

	for (loop = 0; loop < 31 * Data_Acnt; loop++)
		for (loop2 = 0; loop2 < 31 * Data_Bcnt; loop2++)
			if (AAb[loop] == BBb[loop2])
				for (loop3 = 0; loop3 < 31 * Data_Ccnt; loop3++)
					if (AAb[loop] == CCb[loop3])
					{
						theta_temp1[theta_temp1_cnt] = AAb[loop];
						theta_temp1_cnt++;
					}

	theta_temp1_cnt = sort(theta_temp1, theta_temp1_cnt);
	*piangle = p_i / M_PI * 180;

	loop_temp2 = 0;

	for (loop = 0; loop <= theta_temp1_cnt; loop++)
	{
		theta_diff1[loop] = theta_temp1[loop + 1] - theta_temp1[loop];
		if (theta_temp1[loop + 1] - theta_temp1[loop] > 1)
		{
			theta_diff1_idx[loop_temp2] = loop;
			loop_temp2++;
		}
		theta_diff1_idx[loop_temp2] = loop;//마지막까지 참고하기 위함
	}

	temp = 0;
	for (loop = 0; loop <= (int)theta_diff1_idx[0] + 1; loop++)
		temp = temp + theta_temp1[loop];

	temp2 = 0;
	for (loop = (int)theta_diff1_idx[0] + 1; loop <= (int)theta_diff1_idx[1]; loop++)
		temp2 = temp2 + theta_temp1[loop];

	if ((temp / theta_temp1_cnt) > 0)
		*theta = temp / ((int)theta_diff1_idx[0] + 1);

	else
		*theta = temp / ((int)theta_diff1_idx[0] + 1) + 360;

	return 0;
}

double cord3_trans(double *planes_azimuth, double *planes_elevation, double *azimuthes, double *s_piangle_all, int angle_range, double mic_distance)
{
	double osh0[2], osh1[2], osh11[2], osh2[2], osh21[2], osh3[2], osh31[2];
	double osh1_unit1[3], osh1_unit11[3], osh2_unit2[3], osh2_unit21[3], osh3_unit3[3], osh3_unit31[3];
	double a = 0.0, sum = 0.0;
	double temp, temp2[3], temp1[3];
	double tr_out1[3], tr_out2[3], tr_out3[3], tr_out11[3], tr_out21[3], tr_out31[3];
	double dir0[2], dir01[2], dir1[2], dir11[2], dir2[2], dir21[2], dir3[2], dir31[2];
	double x, y, z, t_r1, t_r2, t_r3, t_theta1, t_theta2, t_theta3, t_piangle1, t_piangle2, t_piangle3;
	double tune_theta1 = 0, tune_theta2 = 0, tune_theta3 = 0, tune_pi1 = 0, tune_pi2 = 0, tune_pi3 = 0;
	double dir0_th_array[41], dir0_pi_array[41], dir01_th_array[41], dir01_pi_array[41];
	double dir1_th_array[41], dir1_pi_array[41], dir11_th_array[41], dir11_pi_array[41];
	double dir2_th_array[41], dir2_pi_array[41], dir21_th_array[41], dir21_pi_array[41];
	double dir3_th_array[41], dir3_pi_array[41], dir31_th_array[41], dir31_pi_array[41];
	double dir_th_0p_1p, dir_pi_0p_1p, dir_th_0p_1m, dir_pi_0p_1m, dir_th_0p_2p, dir_pi_0p_2p, dir_th_0p_2m, dir_pi_0p_2m;
	double dir_th_0p_3p, dir_pi_0p_3p, dir_th_0p_3m, dir_pi_0p_3m, dir_th_0m_1p, dir_pi_0m_1p, dir_th_0m_1m, dir_pi_0m_1m;
	double dir_th_0m_2p, dir_pi_0m_2p, dir_th_0m_2m, dir_pi_0m_2m, dir_th_0m_3p, dir_pi_0m_3p, dir_th_0m_3m, dir_pi_0m_3m;
	double dir_th_1p_2p, dir_pi_1p_2p, dir_th_1p_2m, dir_pi_1p_2m, dir_th_1p_3p, dir_pi_1p_3p, dir_th_1p_3m, dir_pi_1p_3m;
	double dir_th_1m_2p, dir_pi_1m_2p, dir_th_1m_2m, dir_pi_1m_2m, dir_th_1m_3p, dir_pi_1m_3p, dir_th_1m_3m, dir_pi_1m_3m;
	double dir_th_2p_3p, dir_pi_2p_3p, dir_th_2p_3m, dir_pi_2p_3m, dir_th_2m_3p, dir_pi_2m_3p, dir_th_2m_3m, dir_pi_2m_3m;
	double th_array[24], pi_array[24];
	double osh_th = 0.0, osh_pi = 0.0;
	int i, k, data_tune, data_len = 0, ret_cnt = 0;

	double y_t[3][3], z_t[3][3], z_tm[3][3], a_axisp[3], a_axism[3];

	//좌표변환 시작
	y_t[0][0] = 0.3333;              y_t[0][1] = 0;                  y_t[0][2] = -0.942809041582063;
	y_t[1][0] = 0;                   y_t[1][1] = 1;                  y_t[1][2] = 0;
	y_t[2][0] = 0.942809041582063;   y_t[2][1] = 0;                  y_t[2][2] = 0.3333;

	z_t[0][0] = -0.5;                z_t[0][1] = 0.866025403784439;  z_t[0][2] = 0;
	z_t[1][0] = -0.866025403784439;  z_t[1][1] = -0.5;               z_t[1][2] = 0;
	z_t[2][0] = 0;                   z_t[2][1] = 0;                  z_t[2][2] = 1;

	z_tm[0][0] = -0.5;               z_tm[0][1] = -0.866025403784439;  z_tm[0][2] = 0;
	z_tm[1][0] = 0.866025403784439;  z_tm[1][1] = -0.5;                z_tm[1][2] = 0;
	z_tm[2][0] = 0;                  z_tm[2][1] = 0;                   z_tm[2][2] = 1;

	osh0[0] = planes_azimuth[0]; 
	osh0[1] = planes_elevation[0];

	osh1[0] = planes_azimuth[1]; 
	osh1[1] = planes_elevation[1];

	osh1_unit1[0] = 100 * (cos(M_PI*osh1[1] / 180)*cos(M_PI*osh1[0] / 180));
	osh1_unit1[1] = 100 * (cos(M_PI*osh1[1] / 180)*sin(M_PI*osh1[0] / 180));
	osh1_unit1[2] = 100 * (sin(M_PI*osh1[1] / 180));

	osh11[0] = osh1[0]; 
	osh11[1] = -1 * osh1[1];
	osh1_unit11[0] = 100 * (cos(M_PI*osh11[1] / 180)*cos(M_PI*osh11[0] / 180));
	osh1_unit11[1] = 100 * (cos(M_PI*osh11[1] / 180)*sin(M_PI*osh11[0] / 180));
	osh1_unit11[2] = 100 * (sin(M_PI*osh11[1] / 180));

	osh2[0] = planes_azimuth[2]; 
	osh2[1] = planes_elevation[2];
	osh2_unit2[0] = 100 * (cos(M_PI*osh2[1] / 180)*cos(M_PI*osh2[0] / 180));
	osh2_unit2[1] = 100 * (cos(M_PI*osh2[1] / 180)*sin(M_PI*osh2[0] / 180));
	osh2_unit2[2] = 100 * (sin(M_PI*osh2[1] / 180));

	osh21[0] = osh2[0]; 
	osh21[1] = -1 * osh2[1];
	osh2_unit21[0] = 100 * (cos(M_PI*osh21[1] / 180)*cos(M_PI*osh21[0] / 180));
	osh2_unit21[1] = 100 * (cos(M_PI*osh21[1] / 180)*sin(M_PI*osh21[0] / 180));
	osh2_unit21[2] = 100 * (sin(M_PI*osh21[1] / 180));

	osh3[0] = planes_azimuth[3];
	osh3[1] = planes_elevation[3];
	osh3_unit3[0] = 100 * (cos(M_PI*osh3[1] / 180)*cos(M_PI*osh3[0] / 180));
	osh3_unit3[1] = 100 * (cos(M_PI*osh3[1] / 180)*sin(M_PI*osh3[0] / 180));
	osh3_unit3[2] = 100 * (sin(M_PI*osh3[1] / 180));

	osh31[0] = osh3[0]; 
	osh31[1] = -1 * osh3[1];
	osh3_unit31[0] = 100 * (cos(M_PI*osh31[1] / 180)*cos(M_PI*osh31[0] / 180));
	osh3_unit31[1] = 100 * (cos(M_PI*osh31[1] / 180)*sin(M_PI*osh31[0] / 180));
	osh3_unit31[2] = 100 * (sin(M_PI*osh31[1] / 180));

	a = mic_distance / sqrt(3.0);
	a_axisp[0] = a / 2; 
	a_axisp[1] = 0;
	a_axisp[2] = 0;
	a_axism[0] = -1 * a / 2;
	a_axism[1] = 0; 
	a_axism[2] = 0;

	// tr_out1/////////////////////////////////////
	temp1[0] = osh1_unit1[0] + a_axisp[0];
	temp1[1] = osh1_unit1[1] + a_axisp[1];
	temp1[2] = osh1_unit1[2] + a_axisp[2];

	for (i = 0; i<3; i++) 
	{
		sum = 0;
		for (k = 0; k<3; k++) 
			sum += y_t[i][k] * temp1[k];
		
		temp2[i] = sum;
	}

	temp1[0] = temp2[0] + a_axism[0];
	temp1[1] = temp2[1] + a_axism[1];
	temp1[2] = temp2[2] + a_axism[2];

	for (i = 0; i<3; i++) 
	{
		sum = 0;
		for (k = 0; k<3; k++) 
			sum += z_t[i][k] * temp1[k];
		
		tr_out1[i] = sum;
	}

	// tr_out2/////////////////////////////////////		
	temp1[0] = osh2_unit2[0] + a_axisp[0];
	temp1[1] = osh2_unit2[1] + a_axisp[1];
	temp1[2] = osh2_unit2[2] + a_axisp[2];

	for (i = 0; i<3; i++)
	{
		sum = 0;
		for (k = 0; k<3; k++) 
			sum += y_t[i][k] * temp1[k];
		
		temp2[i] = sum;
	}

	tr_out2[0] = temp2[0] + a_axism[0];
	tr_out2[1] = temp2[1] + a_axism[1];
	tr_out2[2] = temp2[2] + a_axism[2];

	// tr_out3/////////////////////////////////////
	temp1[0] = osh3_unit3[0] + a_axisp[0];
	temp1[1] = osh3_unit3[1] + a_axisp[1];
	temp1[2] = osh3_unit3[2] + a_axisp[2];

	for (i = 0; i<3; i++) 
	{
		sum = 0;
		for (k = 0; k<3; k++) 
			sum += y_t[i][k] * temp1[k];
		
		temp2[i] = sum;
	}

	temp1[0] = temp2[0] + a_axism[0];
	temp1[1] = temp2[1] + a_axism[1];
	temp1[2] = temp2[2] + a_axism[2];

	for (i = 0; i<3; i++) 
	{
		sum = 0;
		for (k = 0; k<3; k++) 
			sum += z_tm[i][k] * temp1[k];
		
		tr_out3[i] = sum;
	}

	// tr_out11/////////////////////////////////////
	temp1[0] = osh1_unit11[0] + a_axisp[0];
	temp1[1] = osh1_unit11[1] + a_axisp[1];
	temp1[2] = osh1_unit11[2] + a_axisp[2];

	for (i = 0; i<3; i++) 
	{
		sum = 0;
		for (k = 0; k<3; k++) 
			sum += y_t[i][k] * temp1[k];
		
		temp2[i] = sum;
	}

	temp1[0] = temp2[0] + a_axism[0];
	temp1[1] = temp2[1] + a_axism[1];
	temp1[2] = temp2[2] + a_axism[2];

	for (i = 0; i<3; i++) 
	{
		sum = 0;
		for (k = 0; k<3; k++) 
			sum += z_t[i][k] * temp1[k];
		
		tr_out11[i] = sum;
	}

	// tr_out21/////////////////////////////////////		
	temp1[0] = osh2_unit21[0] + a_axisp[0];
	temp1[1] = osh2_unit21[1] + a_axisp[1];
	temp1[2] = osh2_unit21[2] + a_axisp[2];

	for (i = 0; i<3; i++) 
	{
		sum = 0;
		for (k = 0; k<3; k++) 
			sum += y_t[i][k] * temp1[k];
		
		temp2[i] = sum;
	}

	tr_out21[0] = temp2[0] + a_axism[0];
	tr_out21[1] = temp2[1] + a_axism[1];
	tr_out21[2] = temp2[2] + a_axism[2];

	// tr_out31/////////////////////////////////////
	temp1[0] = osh3_unit31[0] + a_axisp[0];
	temp1[1] = osh3_unit31[1] + a_axisp[1];
	temp1[2] = osh3_unit31[2] + a_axisp[2];

	for (i = 0; i<3; i++) 
	{
		sum = 0;
		for (k = 0; k<3; k++) 
			sum += y_t[i][k] * temp1[k];
		
		temp2[i] = sum;
	}

	temp1[0] = temp2[0] + a_axism[0];
	temp1[1] = temp2[1] + a_axism[1];
	temp1[2] = temp2[2] + a_axism[2];

	for (i = 0; i < 3; i++) 
	{
		sum = 0;
		for (k = 0; k < 3; k++) 
			sum += z_tm[i][k] * temp1[k];
		
		tr_out31[i] = sum;
	}

	dir0[0] = osh0[0]; 
	dir0[1] = osh0[1];
	dir01[0] = osh0[0]; 
	dir01[1] = -osh0[1];

	x = tr_out1[0];           
	y = tr_out1[1];            
	z = tr_out1[2];

	t_r1 = sqrt(x*x + y*y + z*z);
	t_theta1 = (acos(x / sqrt(x*x + y*y)) / (2 * M_PI)) * 360;

	if (y < 0) 
		t_theta1 = 360 - t_theta1;
	
	t_piangle1 = 90 - (atan(sqrt(x*x + y*y) / z) / M_PI) * 180;

	if (t_piangle1 > 90) 
		t_piangle1 = 180 - t_piangle1;
	
	dir1[0] = round(t_theta1);
	dir1[1] = round(t_piangle1);

	x = tr_out11[0];         
	y = tr_out11[1];            
	z = tr_out11[2];

	t_r1 = sqrt(x*x + y*y + z*z);
	t_theta1 = (acos(x / sqrt(x*x + y*y)) / (2 * M_PI)) * 360;
	if (y < 0) 
		t_theta1 = 360 - t_theta1;
	
	t_piangle1 = 90 - (atan(sqrt(x*x + y*y) / z) / M_PI) * 180;
	if (t_piangle1 > 90) 
		t_piangle1 = 180 - t_piangle1;
	
	dir11[0] = round(t_theta1);
	dir11[1] = round(t_piangle1);

	x = tr_out2[0];            
	y = tr_out2[1];       
	z = tr_out2[2];

	t_r2 = sqrt(x*x + y*y + z*z);
	t_theta2 = (acos(x / sqrt(x*x + y*y)) / (2 * M_PI)) * 360;
	if (y < 0) 
		t_theta2 = 360 - t_theta2;
	
	t_piangle2 = 90 - (atan(sqrt(x*x + y*y) / z) / M_PI) * 180;
	if (t_piangle2 > 90) 
		t_piangle2 = 180 - t_piangle2;
	
	dir2[0] = round(t_theta2);
	dir2[1] = round(t_piangle2);

	x = tr_out21[0];          
	y = tr_out21[1];        
	z = tr_out21[2];

	t_r2 = sqrt(x*x + y*y + z*z);
	t_theta2 = (acos(x / sqrt(x*x + y*y)) / (2 * M_PI)) * 360;
	if (y < 0) 
		t_theta2 = 360 - t_theta2;
	
	t_piangle2 = 90 - (atan(sqrt(x*x + y*y) / z) / M_PI) * 180;
	if (t_piangle2 > 90) 
		t_piangle2 = 180 - t_piangle2;
	
	dir21[0] = round(t_theta2);
	dir21[1] = round(t_piangle2);

	x = tr_out3[0];           
	y = tr_out3[1];      
	z = tr_out3[2];

	t_r3 = sqrt(x*x + y*y + z*z);
	t_theta3 = (acos(x / sqrt(x*x + y*y)) / (2 * M_PI)) * 360;

	if (y < 0) 
		t_theta3 = 360 - t_theta3;
	
	t_piangle3 = 90 - (atan(sqrt(x*x + y*y) / z) / M_PI) * 180;
	if (t_piangle3 > 90) 
		t_piangle3 = 180 - t_piangle3;
	
	dir3[0] = round(t_theta3);
	dir3[1] = round(t_piangle3);

	x = tr_out31[0];           
	y = tr_out31[1];        
	z = tr_out31[2];

	t_r3 = sqrt(x*x + y*y + z*z);
	t_theta3 = (acos(x / sqrt(x*x + y*y)) / (2 * M_PI)) * 360;

	if (y < 0)  // % 뭔가틀린건데 모르겠음.. y부호가 반대임
		t_theta3 = 360 - t_theta3;
	
	t_piangle3 = 90 - (atan(sqrt(x*x + y*y) / z) / M_PI) * 180;
	if (t_piangle3 > 90) 
		t_piangle3 = 180 - t_piangle3;
	
	dir31[0] = round(t_theta3);
	dir31[1] = round(t_piangle3);

	dir1[0] = dir1[0] + tune_theta1;
	dir1[1] = dir1[1] + tune_pi1;

	dir11[0] = dir11[0] + tune_theta1;
	dir11[1] = dir11[1] + tune_pi1;

	dir2[0] = dir2[0] + tune_theta2;
	dir2[1] = dir2[1] + tune_pi2;

	dir21[0] = dir21[0] + tune_theta2;
	dir21[1] = dir21[1] + tune_pi2;

	dir3[0] = dir3[0] + tune_theta3;
	dir3[1] = dir3[1] + tune_pi3;

	dir31[0] = dir31[0] + tune_theta3;
	dir31[1] = dir31[1] + tune_pi3;
	
	//좌표변환 완료 

	data_tune = angle_range;
	for (i = -1 * data_tune; i <= data_tune; i++)
	{
		dir0_th_array[i + data_tune] = round(dir0[0] + i);
		dir0_pi_array[i + data_tune] = round(dir0[1] + i);

		dir01_th_array[i + data_tune] = round(dir01[0] + i);
		dir01_pi_array[i + data_tune] = round(dir01[1] + i);

		dir1_th_array[i + data_tune] = round(dir1[0] + i);
		dir1_pi_array[i + data_tune] = round(dir1[1] + i);

		dir11_th_array[i + data_tune] = round(dir11[0] + i);
		dir11_pi_array[i + data_tune] = round(dir11[1] + i);

		dir2_th_array[i + data_tune] = round(dir2[0] + i);
		dir2_pi_array[i + data_tune] = round(dir2[1] + i);

		dir21_th_array[i + data_tune] = round(dir21[0] + i);
		dir21_pi_array[i + data_tune] = round(dir21[1] + i);

		dir3_th_array[i + data_tune] = round(dir3[0] + i);
		dir3_pi_array[i + data_tune] = round(dir3[1] + i);

		dir31_th_array[i + data_tune] = round(dir31[0] + i);
		dir31_pi_array[i + data_tune] = round(dir31[1] + i);
	}

	data_len = (data_tune * 2) + 1;
	temp = sortd(dir0_th_array, data_len);
	temp = sortd(dir0_pi_array, data_len);

	temp = sortd(dir01_th_array, data_len);
	temp = sortd(dir01_pi_array, data_len);

	temp = sortd(dir1_th_array, data_len);
	temp = sortd(dir1_pi_array, data_len);

	temp = sortd(dir11_th_array, data_len);
	temp = sortd(dir11_pi_array, data_len);

	temp = sortd(dir2_th_array, data_len);
	temp = sortd(dir2_pi_array, data_len);

	temp = sortd(dir21_th_array, data_len);
	temp = sortd(dir21_pi_array, data_len);

	temp = sortd(dir3_th_array, data_len);
	temp = sortd(dir3_pi_array, data_len);

	temp = sortd(dir31_th_array, data_len);
	temp = sortd(dir31_pi_array, data_len);

	// 교차값이 없으면 리턴값 NO_INTERSECT_VALUE			
	dir_th_0p_1p = intsec_mean(dir0_th_array, dir1_th_array, data_len);
	dir_pi_0p_1p = intsec_mean(dir0_pi_array, dir1_pi_array, data_len);
	if (dir_th_0p_1p != NO_INTERSECT_VALUE && dir_pi_0p_1p != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_0p_1p;
		pi_array[ret_cnt] = dir_pi_0p_1p;
		ret_cnt++;
	}

	dir_th_0p_1m = intsec_mean(dir0_th_array, dir11_th_array, data_len);
	dir_pi_0p_1m = intsec_mean(dir0_pi_array, dir11_pi_array, data_len);
	if (dir_th_0p_1m != NO_INTERSECT_VALUE && dir_pi_0p_1m != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_0p_1m;
		pi_array[ret_cnt] = dir_pi_0p_1m;
		ret_cnt++;
	}

	dir_th_0p_2p = intsec_mean(dir0_th_array, dir2_th_array, data_len);
	dir_pi_0p_2p = intsec_mean(dir0_pi_array, dir2_pi_array, data_len);
	if (dir_th_0p_2p != NO_INTERSECT_VALUE && dir_pi_0p_2p != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_0p_2p;
		pi_array[ret_cnt] = dir_pi_0p_2p;
		ret_cnt++;
	}

	dir_th_0p_2m = intsec_mean(dir0_th_array, dir21_th_array, data_len);
	dir_pi_0p_2m = intsec_mean(dir0_pi_array, dir21_pi_array, data_len);
	if (dir_th_0p_2m != NO_INTERSECT_VALUE && dir_pi_0p_2m != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_0p_2m;
		pi_array[ret_cnt] = dir_pi_0p_2m;
		ret_cnt++;
	}

	dir_th_0p_3p = intsec_mean(dir0_th_array, dir3_th_array, data_len);
	dir_pi_0p_3p = intsec_mean(dir0_pi_array, dir3_pi_array, data_len);
	if (dir_th_0p_3p != NO_INTERSECT_VALUE && dir_pi_0p_3p != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_0p_3p;
		pi_array[ret_cnt] = dir_pi_0p_3p;
		ret_cnt++;
	}

	dir_th_0p_3m = intsec_mean(dir0_th_array, dir31_th_array, data_len);
	dir_pi_0p_3m = intsec_mean(dir0_pi_array, dir31_pi_array, data_len);
	if (dir_th_0p_3m != NO_INTERSECT_VALUE && dir_pi_0p_3m != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_0p_3m;
		pi_array[ret_cnt] = dir_pi_0p_3m;
		ret_cnt++;
	}
	/////////////////////////////////////////

	dir_th_0m_1p = intsec_mean(dir01_th_array, dir1_th_array, data_len);
	dir_pi_0m_1p = intsec_mean(dir01_pi_array, dir1_pi_array, data_len);
	if (dir_th_0m_1p != NO_INTERSECT_VALUE && dir_pi_0m_1p != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_0m_1p;
		pi_array[ret_cnt] = dir_pi_0m_1p;
		ret_cnt++;
	}

	dir_th_0m_1m = intsec_mean(dir01_th_array, dir11_th_array, data_len);
	dir_pi_0m_1m = intsec_mean(dir01_pi_array, dir11_pi_array, data_len);
	if (dir_th_0m_1m != NO_INTERSECT_VALUE && dir_pi_0m_1m != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_0m_1m;                  
		pi_array[ret_cnt] = dir_pi_0m_1m;
		ret_cnt++;
	}

	dir_th_0m_2p = intsec_mean(dir01_th_array, dir2_th_array, data_len);
	dir_pi_0m_2p = intsec_mean(dir01_pi_array, dir2_pi_array, data_len);
	if (dir_th_0m_2p != NO_INTERSECT_VALUE && dir_pi_0m_2p != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_0m_2p;               
		pi_array[ret_cnt] = dir_pi_0m_2p;
		ret_cnt++;
	}

	dir_th_0m_2m = intsec_mean(dir01_th_array, dir21_th_array, data_len);
	dir_pi_0m_2m = intsec_mean(dir01_pi_array, dir21_pi_array, data_len);
	if (dir_th_0m_2m != NO_INTERSECT_VALUE && dir_pi_0m_2m != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_0m_2m;                  
		pi_array[ret_cnt] = dir_pi_0m_2m;
		ret_cnt++;
	}

	dir_th_0m_3p = intsec_mean(dir01_th_array, dir3_th_array, data_len);
	dir_pi_0m_3p = intsec_mean(dir01_pi_array, dir3_pi_array, data_len);
	if (dir_th_0m_3p != NO_INTERSECT_VALUE && dir_pi_0m_3p != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_0m_3p;                
		pi_array[ret_cnt] = dir_pi_0m_3p;
		ret_cnt++;
	}

	dir_th_0m_3m = intsec_mean(dir01_th_array, dir31_th_array, data_len);
	dir_pi_0m_3m = intsec_mean(dir01_pi_array, dir31_pi_array, data_len);
	if (dir_th_0m_3m != NO_INTERSECT_VALUE && dir_pi_0m_3m != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_0m_3m;                   
		pi_array[ret_cnt] = dir_pi_0m_3m;
		ret_cnt++;
	}
	///////////////////////////////////////

	dir_th_1p_2p = intsec_mean(dir1_th_array, dir2_th_array, data_len);
	dir_pi_1p_2p = intsec_mean(dir1_pi_array, dir2_pi_array, data_len);
	if (dir_th_1p_2p != NO_INTERSECT_VALUE && dir_pi_1p_2p != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_1p_2p;            
		pi_array[ret_cnt] = dir_pi_1p_2p;
		ret_cnt++;
	}

	dir_th_1p_2m = intsec_mean(dir1_th_array, dir21_th_array, data_len);
	dir_pi_1p_2m = intsec_mean(dir1_pi_array, dir21_pi_array, data_len);
	if (dir_th_1p_2m != NO_INTERSECT_VALUE && dir_pi_1p_2m != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_1p_2m;               
		pi_array[ret_cnt] = dir_pi_1p_2m;
		ret_cnt++;
	}
	dir_th_1p_3p = intsec_mean(dir1_th_array, dir3_th_array, data_len);
	dir_pi_1p_3p = intsec_mean(dir1_pi_array, dir3_pi_array, data_len);
	if (dir_th_1p_3p != NO_INTERSECT_VALUE && dir_pi_1p_3p != NO_INTERSECT_VALUE) 
	{
		th_array[ret_cnt] = dir_th_1p_3p;                  
		pi_array[ret_cnt] = dir_pi_1p_3p;
		ret_cnt++;
	}

	dir_th_1p_3m = intsec_mean(dir1_th_array, dir31_th_array, data_len);
	dir_pi_1p_3m = intsec_mean(dir1_pi_array, dir31_pi_array, data_len);
	if (dir_th_1p_3m != NO_INTERSECT_VALUE && dir_pi_1p_3m != NO_INTERSECT_VALUE) 
	{
		th_array[ret_cnt] = dir_th_1p_3m;                  
		pi_array[ret_cnt] = dir_pi_1p_3m;
		ret_cnt++;
	}

	dir_th_1m_2p = intsec_mean(dir11_th_array, dir2_th_array, data_len);
	dir_pi_1m_2p = intsec_mean(dir11_pi_array, dir2_pi_array, data_len);
	if (dir_th_1m_2p != NO_INTERSECT_VALUE && dir_pi_1m_2p != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_1m_2p;                   
		pi_array[ret_cnt] = dir_pi_1m_2p;
		ret_cnt++;
	}

	dir_th_1m_2m = intsec_mean(dir11_th_array, dir21_th_array, data_len);
	dir_pi_1m_2m = intsec_mean(dir11_pi_array, dir21_pi_array, data_len);
	if (dir_th_1m_2m != NO_INTERSECT_VALUE && dir_pi_1m_2m != NO_INTERSECT_VALUE) 
	{
		th_array[ret_cnt] = dir_th_1m_2m;                 
		pi_array[ret_cnt] = dir_pi_1m_2m;
		ret_cnt++;
	}

	dir_th_1m_3p = intsec_mean(dir11_th_array, dir3_th_array, data_len);
	dir_pi_1m_3p = intsec_mean(dir11_pi_array, dir3_pi_array, data_len);

	if (dir_th_1m_3p != NO_INTERSECT_VALUE && dir_pi_1m_3p != NO_INTERSECT_VALUE) 
	{
		th_array[ret_cnt] = dir_th_1m_3p;            
		pi_array[ret_cnt] = dir_pi_1m_3p;
		ret_cnt++;
	}

	dir_th_1m_3m = intsec_mean(dir11_th_array, dir31_th_array, data_len);
	dir_pi_1m_3m = intsec_mean(dir11_pi_array, dir31_pi_array, data_len);
	if (dir_th_1m_3m != NO_INTERSECT_VALUE && dir_pi_1m_3m != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_1m_3m;               
		pi_array[ret_cnt] = dir_pi_1m_3m;
		ret_cnt++;
	}

	dir_th_2p_3p = intsec_mean(dir2_th_array, dir3_th_array, data_len);
	dir_pi_2p_3p = intsec_mean(dir2_pi_array, dir3_pi_array, data_len);
	if (dir_th_2p_3p != NO_INTERSECT_VALUE && dir_pi_2p_3p != NO_INTERSECT_VALUE) 
	{
		th_array[ret_cnt] = dir_th_2p_3p;                
		pi_array[ret_cnt] = dir_pi_2p_3p;
		ret_cnt++;
	}

	dir_th_2p_3m = intsec_mean(dir2_th_array, dir31_th_array, data_len);
	dir_pi_2p_3m = intsec_mean(dir2_pi_array, dir31_pi_array, data_len);
	if (dir_th_2p_3m != NO_INTERSECT_VALUE && dir_pi_2p_3m != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_2p_3m;                
		pi_array[ret_cnt] = dir_pi_2p_3m;
		ret_cnt++;
	}

	dir_th_2m_3p = intsec_mean(dir21_th_array, dir3_th_array, data_len);
	dir_pi_2m_3p = intsec_mean(dir21_pi_array, dir3_pi_array, data_len);
	if (dir_th_2m_3p != NO_INTERSECT_VALUE && dir_pi_2m_3p != NO_INTERSECT_VALUE)
	{
		th_array[ret_cnt] = dir_th_2m_3p;                  
		pi_array[ret_cnt] = dir_pi_2m_3p;
		ret_cnt++;
	}

	dir_th_2m_3m = intsec_mean(dir21_th_array, dir31_th_array, data_len);
	dir_pi_2m_3m = intsec_mean(dir21_pi_array, dir31_pi_array, data_len);
	if (dir_th_2m_3m != NO_INTERSECT_VALUE && dir_pi_2m_3m != NO_INTERSECT_VALUE) 
	{
		th_array[ret_cnt] = dir_th_2m_3m;                  
		pi_array[ret_cnt] = dir_pi_2m_3m;
		ret_cnt++;
	}

	// 4개의 각들 평균 냄
	for (i = 0; i < ret_cnt; i++) 
	{
		osh_th = osh_th + th_array[i];					
		osh_pi = osh_pi + pi_array[i];
	}

	planes_azimuth[4] = osh_th / (double) ret_cnt;
	planes_elevation[4] = osh_pi / (double)ret_cnt;

	if (planes_azimuth[4] == 360)
		planes_azimuth[4] = 0;

	azimuthes[0] = dir0[0];
	azimuthes[1] = dir01[0];
	azimuthes[2] = dir1[0];
	azimuthes[3] = dir11[0];
	azimuthes[4] = dir2[0];
	azimuthes[5] = dir21[0];
	azimuthes[6] = dir3[0];
	azimuthes[7] = dir31[0];

	//printf("%lf %lf %lf %lf %lf %lf %lf %lf\n", azimuthes[0], azimuthes[1], azimuthes[2], azimuthes[3], azimuthes[4],
	//	azimuthes[5], azimuthes[6], azimuthes[7]);

	s_piangle_all[0] = dir0[1];
	s_piangle_all[1] = dir01[1];
	s_piangle_all[2] = dir1[1];
	s_piangle_all[3] = dir11[1];
	s_piangle_all[4] = dir2[1];
	s_piangle_all[5] = dir21[1];
	s_piangle_all[6] = dir3[1];
	s_piangle_all[7] = dir31[1];

	return temp;
}
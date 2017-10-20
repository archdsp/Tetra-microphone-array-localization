#pragma once

#define SIGNAL_BUFF 30000 /*gs_buff_s �� 3�� �޾Ƽ� ó����*/

// gs_buff_s �� 3�� �޾Ƽ� ó����
#define FFT_BUFF 32768 

typedef struct
{
	unsigned int length;       /* size of filter */
	double *coef;               /* pointer to coefficients of filter */
} FILTER;


void cdft(int n, int isgn, double *a);
void fir_filter_array(double *in, double *out, int in_len, FILTER *fir);
double CPSP_FILT_normal(double *Buff0, double *Buff1, int Len, double max_d);

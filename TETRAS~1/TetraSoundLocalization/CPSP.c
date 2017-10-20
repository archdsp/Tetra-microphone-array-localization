#pragma once
#include <math.h>
#include <stdio.h>
#include "CPSP.h"

// ILD(Interaural level difference) 를 위한 HPF
const double hpf[257] = { 0.0001533,0.0001506,0.0001483,0.0001463,0.0001444,0.0001428,0.0001413,0.0001399,0.0001385,0.0001371,					0.0001355,0.0001338,0.0001319,0.0001297,0.0001271,0.0001240,0.0001203,0.0001160,0.0001110,0.0001052,
0.0000984,0.0000907,0.0000818,0.0000718,0.0000604,0.0000477,0.0000334,0.0000176,0.0000000,-0.0000194,				-0.0000406,-0.0000639,-0.0000892,-0.0001168,-0.0001466,-0.0001788,-0.0002135,-0.0002508,-0.0002906,-0.0003333,
-0.0003787,-0.0004270,-0.0004782,-0.0005325,-0.0005899,-0.0006503,-0.0007140,-0.0007809,-0.0008511,-0.0009246,		-0.0010014,-0.0010816,-0.0011652,-0.0012522,-0.0013426,-0.0014364,-0.0015335,-0.0016341,-0.0017380,-0.0018452,
-0.0019557,-0.0020695,-0.0021865,-0.0023066,-0.0024298,-0.0025559,-0.0026850,-0.0028169,-0.0029516,-0.0030889,		-0.0032288,-0.0033710,-0.0035156,-0.0036623,-0.0038111,-0.0039617,-0.0041141,-0.0042682,-0.0044236,-0.0045804,
-0.0047382,-0.0048971,-0.0050567,-0.0052168,-0.0053774,-0.0055382,-0.0056991,-0.0058598,-0.0060201,-0.0061799,		-0.0063389,-0.0064970,-0.0066539,-0.0068095,-0.0069635,-0.0071157,-0.0072660,-0.0074141,-0.0075599,-0.0077031,
-0.0078435,-0.0079810,-0.0081154,-0.0082464,-0.0083739,-0.0084978,-0.0086177,-0.0087337,-0.0088454,-0.0089528,		-0.0090556,-0.0091538,-0.0092472,-0.0093356,-0.0094189,-0.0094971,-0.0095699,-0.0096373,-0.0096992,-0.0097555,
-0.0098061,-0.0098509,-0.0098899,-0.0099229,-0.0099501,-0.0099712,-0.0099863,-0.0099954,0.9898469,-0.0099954,		-0.0099863,-0.0099712,-0.0099501,-0.0099229,-0.0098899,-0.0098509,-0.0098061,-0.0097555,-0.0096992,-0.0096373,
-0.0095699,-0.0094971,-0.0094189,-0.0093356,-0.0092472,-0.0091538,-0.0090556,-0.0089528,-0.0088454,-0.0087337,		-0.0086177,-0.0084978,-0.0083739,-0.0082464,-0.0081154,-0.0079810,-0.0078435,-0.0077031,-0.0075599,-0.0074141,
-0.0072660,-0.0071157,-0.0069635,-0.0068095,-0.0066539,-0.0064970,-0.0063389,-0.0061799,-0.0060201,-0.0058598,		-0.0056991,-0.0055382,-0.0053774,-0.0052168,-0.0050567,-0.0048971,-0.0047382,-0.0045804,-0.0044236,-0.0042682,
-0.0041141,-0.0039617,-0.0038111,-0.0036623,-0.0035156,-0.0033710,-0.0032288,-0.0030889,-0.0029516,-0.0028169,		-0.0026850,-0.0025559,-0.0024298,-0.0023066,-0.0021865,-0.0020695,-0.0019557,-0.0018452,-0.0017380,-0.0016341,
-0.0015335,-0.0014364,-0.0013426,-0.0012522,-0.0011652,-0.0010816,-0.0010014,-0.0009246,-0.0008511,-0.0007809,		-0.0007140,-0.0006503,-0.0005899,-0.0005325,-0.0004782,-0.0004270,-0.0003787,-0.0003333,-0.0002906,-0.0002508,
-0.0002135,-0.0001788,-0.0001466,-0.0001168,-0.0000892,-0.0000639,-0.0000406,-0.0000194,0.0000000,0.0000176,		0.0000334,0.0000477,0.0000604,0.0000718,0.0000818,0.0000907,0.0000984,0.0001052,0.0001110,0.0001160,
0.0001203,0.0001240,0.0001271,0.0001297,0.0001319,0.0001338,0.0001355,0.0001371,0.0001385,0.0001399,				0.0001413,0.0001428,0.0001444,0.0001463,0.0001483,0.0001506,0.0001533 };

double	m_pExecBuffer_t[2][FFT_BUFF * 2];
double	m_pSizeBuffer_1[2][FFT_BUFF * 2];
double	m_pResult[FFT_BUFF * 2];
double	m_pDelay[FFT_BUFF * 2];
double	m_pRLxy[FFT_BUFF * 2];

double CPSP_FILT_normal(double *Buff0, double *Buff1, int Len, double max_d) 
{
	int i;
	double maxindex = 0, TheMaxPos = 0;
	double max = -10000001.0;

	// 데이터에서 LEFT and RIGHT 신호 분리 8단계 1/16
	for (i = 0; i < SIGNAL_BUFF; i++)
	{
		m_pExecBuffer_t[0][i * 2] = Buff0[i];
		m_pExecBuffer_t[0][i * 2 + 1] = 0;
		m_pExecBuffer_t[1][i * 2] = Buff1[i];
		m_pExecBuffer_t[1][i * 2 + 1] = 0;
	}

	// 패딩
	for (; i < FFT_BUFF; i++)
	{
		m_pExecBuffer_t[0][i * 2] = 0;
		m_pExecBuffer_t[0][i * 2 + 1] = 0;
		m_pExecBuffer_t[1][i * 2] = 0;
		m_pExecBuffer_t[1][i * 2 + 1] = 0;
	}

	//    F F T
	cdft(Len * 2, -1, m_pExecBuffer_t[1]);
	cdft(Len * 2, -1, m_pExecBuffer_t[0]);

	//  CPSP 계산 과정.. 1/8
	for (i = 0; i < Len; i++)
	{
		m_pSizeBuffer_1[0][i] = // Auto Power spectrum density를위한 과정
			((m_pExecBuffer_t[0][i * 2] * m_pExecBuffer_t[0][i * 2]
				+ (m_pExecBuffer_t[0][i * 2 + 1]) * (m_pExecBuffer_t[0][i * 2 + 1])));

		m_pSizeBuffer_1[1][i] = // Auto Power spectrum density를위한 과정
			((m_pExecBuffer_t[1][i * 2] * m_pExecBuffer_t[1][i * 2]
				+ (m_pExecBuffer_t[1][i * 2 + 1]) * (m_pExecBuffer_t[1][i * 2 + 1])));

		// Cross Power spectrum density를위한 과정
		m_pRLxy[i * 2] = (m_pExecBuffer_t[0][i * 2] * m_pExecBuffer_t[1][i * 2]
			+ m_pExecBuffer_t[0][i * 2 + 1] * m_pExecBuffer_t[1][i * 2 + 1]);

		// Cross Power spectrum density를위한 과정
		m_pRLxy[i * 2 + 1] = (-1 * m_pExecBuffer_t[0][i * 2 + 1] * m_pExecBuffer_t[1][i * 2]
			+ m_pExecBuffer_t[0][i * 2] * (m_pExecBuffer_t[1][i * 2 + 1]));

		m_pResult[i * 2] = m_pRLxy[i * 2] / sqrt(m_pRLxy[i * 2] * m_pRLxy[i * 2]
			+ m_pRLxy[i * 2 + 1] * m_pRLxy[i * 2 + 1]);

		m_pResult[i * 2 + 1] = m_pRLxy[i * 2 + 1] / sqrt(m_pRLxy[i * 2] * m_pRLxy[i * 2]
			+ m_pRLxy[i * 2 + 1] * m_pRLxy[i * 2 + 1]);

	}

	//  CPSP  에 대한  Inverse Fourier Transform
	cdft(Len * 2, 1, m_pResult);

	for (i = 0; i < Len * 2; i++)
		m_pResult[i] = m_pResult[i] / Len;

	//	FFTSHIFT 과정...
	for (i = 0; i < Len; i++)
	{
		if (i < Len / 2)
			m_pDelay[i * 2] = m_pResult[(i + Len / 2) * 2];

		if (i >= Len / 2)
			m_pDelay[i * 2] = m_pResult[(i - Len / 2) * 2];
	}

	//	실험 결과 저장...
	for (i = (Len / 2) - max_d; i < (Len / 2) + max_d; i++)
		if (m_pDelay[2 * i] > max)
		{
			maxindex = i;
			max = m_pDelay[2 * i];
		}

	//Normalization 
	for (i = 0; i < Len; i++)
		m_pDelay[2 * i] = m_pDelay[2 * i] / max;

	TheMaxPos = (maxindex - Len / 2);

	return TheMaxPos;
}

//  FAST FOURIER TRANSFORM ALGORITHM SUB ROUTINE 여기 부터 끝까지..
/*
Fast Fourier/Cosine/Sine Transform
dimension   :one
data length :power of 2
decimation  :frequency
radix       :4, 2
data        :inplace
table       :not use

functions
cdft: Complex Discrete Fourier Transform
rdft: Real Discrete Fourier Transform
ddct: Discrete Cosine Transform
ddst: Discrete Sine Transform
dfct: Cosine Transform of RDFT (Real Symmetric DFT)
dfst: Sine Transform of RDFT (Real Anti-symmetric DFT)

function prototypes
void cdft(int, int, double *);
void rdft(int, int, double *);
void ddct(int, int, double *);
void ddst(int, int, double *);
void dfct(int, double *);
void dfst(int, double *);


-------- Complex DFT (Discrete Fourier Transform) --------
[definition]
<case1>
X[k] = sum_j=0^n-1 x[j]*exp(2*pi*i*j*k/n), 0<=k<n
<case2>
X[k] = sum_j=0^n-1 x[j]*exp(-2*pi*i*j*k/n), 0<=k<n
(notes: sum_j=0^n-1 is a summation from j=0 to n-1)

[usage]
<case1>
cdft(2*n, 1, a);
<case2>
cdft(2*n, -1, a);
[parameters]
2*n            :data length (int)
n >= 1, n = power of 2
a[0...2*n-1]   :input/output data (double *)
input data
a[2*j] = Re(x[j]),
a[2*j+1] = Im(x[j]), 0<=j<n
output data
a[2*k] = Re(X[k]),
a[2*k+1] = Im(X[k]), 0<=k<n
[remark]
Inverse of
cdft(2*n, -1, a);
is
cdft(2*n, 1, a);
for (j = 0; j <= 2 * n - 1; j++) {
a[j] *= 1.0 / n;
}
.


-------- Real DFT / Inverse of Real DFT --------
[definition]
<case1> RDFT
R[k] = sum_j=0^n-1 a[j]*cos(2*pi*j*k/n), 0<=k<=n/2
I[k] = sum_j=0^n-1 a[j]*sin(2*pi*j*k/n), 0<k<n/2
<case2> IRDFT (excluding scale)
a[k] = (R[0] + R[n/2]*cos(pi*k))/2 +
sum_j=1^n/2-1 R[j]*cos(2*pi*j*k/n) +
sum_j=1^n/2-1 I[j]*sin(2*pi*j*k/n), 0<=k<n
[usage]
<case1>
rdft(n, 1, a);
<case2>
rdft(n, -1, a);
[parameters]
n              :data length (int)
n >= 2, n = power of 2
a[0...n-1]     :input/output data (double *)
<case1>
output data
a[2*k] = R[k], 0<=k<n/2
a[2*k+1] = I[k], 0<k<n/2
a[1] = R[n/2]
<case2>
input data
a[2*j] = R[j], 0<=j<n/2
a[2*j+1] = I[j], 0<j<n/2
a[1] = R[n/2]
[remark]
Inverse of
rdft(n, 1, a);
is
rdft(n, -1, a);
for (j = 0; j <= n - 1; j++) {
a[j] *= 2.0 / n;
}
.


-------- DCT (Discrete Cosine Transform) / Inverse of DCT --------
[definition]
<case1> IDCT (excluding scale)
C[k] = sum_j=0^n-1 a[j]*cos(pi*j*(k+1/2)/n), 0<=k<n
<case2> DCT
C[k] = sum_j=0^n-1 a[j]*cos(pi*(j+1/2)*k/n), 0<=k<n
[usage]
<case1>
ddct(n, 1, a);
<case2>
ddct(n, -1, a);
[parameters]
n              :data length (int)
n >= 2, n = power of 2
a[0...n-1]     :input/output data (double *)
output data
a[k] = C[k], 0<=k<n
[remark]
Inverse of
ddct(n, -1, a);
is
a[0] *= 0.5;
ddct(n, 1, a);
for (j = 0; j <= n - 1; j++) {
a[j] *= 2.0 / n;
}
.


-------- DST (Discrete Sine Transform) / Inverse of DST --------
[definition]
<case1> IDST (excluding scale)
S[k] = sum_j=1^n A[j]*sin(pi*j*(k+1/2)/n), 0<=k<n
<case2> DST
S[k] = sum_j=0^n-1 a[j]*sin(pi*(j+1/2)*k/n), 0<k<=n
[usage]
<case1>
ddst(n, 1, a);
<case2>
ddst(n, -1, a);
[parameters]
n              :data length (int)
n >= 2, n = power of 2
a[0...n-1]     :input/output data (double *)
<case1>
input data
a[j] = A[j], 0<j<n
a[0] = A[n]
output data
a[k] = S[k], 0<=k<n
<case2>
output data
a[k] = S[k], 0<k<n
a[0] = S[n]
[remark]
Inverse of
ddst(n, -1, a);
is
a[0] *= 0.5;
ddst(n, 1, a);
for (j = 0; j <= n - 1; j++) {
a[j] *= 2.0 / n;
}
.


-------- Cosine Transform of RDFT (Real Symmetric DFT) --------
[definition]
C[k] = sum_j=0^n a[j]*cos(pi*j*k/n), 0<=k<=n
[usage]
dfct(n, a);
[parameters]
n              :data length - 1 (int)
n >= 2, n = power of 2
a[0...n]       :input/output data (double *)
output data
a[k] = C[k], 0<=k<=n
[remark]
Inverse of
a[0] *= 0.5;
a[n] *= 0.5;
dfct(n, a);
is
a[0] *= 0.5;
a[n] *= 0.5;
dfct(n, a);
for (j = 0; j <= n; j++) {
a[j] *= 2.0 / n;
}
.


-------- Sine Transform of RDFT (Real Anti-symmetric DFT) --------
[definition]
S[k] = sum_j=1^n-1 a[j]*sin(pi*j*k/n), 0<k<n
[usage]
dfst(n, a);
[parameters]
n              :data length + 1 (int)
n >= 2, n = power of 2
a[0...n-1]     :input/output data (double *)
output data
a[k] = S[k], 0<k<n
(a[0] is used for work area)
[remark]
Inverse of
dfst(n, a);
is
dfst(n, a);
for (j = 1; j <= n - 1; j++) {
a[j] *= 2.0 / n;
}
.
*/


void cdft(int n, int isgn, double *a)
{
	void bitrv2(int n, double *a);
	void bitrv2conj(int n, double *a);
	void cftfsub(int n, double *a);
	void cftbsub(int n, double *a);

	if (n > 4) {
		if (isgn >= 0) {
			bitrv2(n, a);
			cftfsub(n, a);
		}
		else {
			bitrv2conj(n, a);
			cftbsub(n, a);
		}
	}
	else if (n == 4) {
		cftfsub(n, a);
	}
}


void rdft(int n, int isgn, double *a)
{
	void bitrv2(int n, double *a);
	void cftfsub(int n, double *a);
	void cftbsub(int n, double *a);
	void rftfsub(int n, double *a);
	void rftbsub(int n, double *a);
	double xi;

	if (isgn >= 0) {
		if (n > 4) {
			bitrv2(n, a);
			cftfsub(n, a);
			rftfsub(n, a);
		}
		else if (n == 4) {
			cftfsub(n, a);
		}
		xi = a[0] - a[1];
		a[0] += a[1];
		a[1] = xi;
	}
	else {
		a[1] = 0.5 * (a[0] - a[1]);
		a[0] -= a[1];
		if (n > 4) {
			rftbsub(n, a);
			bitrv2(n, a);
			cftbsub(n, a);
		}
		else if (n == 4) {
			cftfsub(n, a);
		}
	}
}


void ddct(int n, int isgn, double *a)
{
	void bitrv2(int n, double *a);
	void cftfsub(int n, double *a);
	void cftbsub(int n, double *a);
	void rftfsub(int n, double *a);
	void rftbsub(int n, double *a);
	void dctsub(int n, double *a);
	void dctsub4(int n, double *a);
	int j;
	double xr;

	if (isgn < 0) {
		xr = a[n - 1];
		for (j = n - 2; j >= 2; j -= 2) {
			a[j + 1] = a[j] - a[j - 1];
			a[j] += a[j - 1];
		}
		a[1] = a[0] - xr;
		a[0] += xr;
		if (n > 4) {
			rftbsub(n, a);
			bitrv2(n, a);
			cftbsub(n, a);
		}
		else if (n == 4) {
			cftfsub(n, a);
		}
	}
	if (n > 4) {
		dctsub(n, a);
	}
	else {
		dctsub4(n, a);
	}
	if (isgn >= 0) {
		if (n > 4) {
			bitrv2(n, a);
			cftfsub(n, a);
			rftfsub(n, a);
		}
		else if (n == 4) {
			cftfsub(n, a);
		}
		xr = a[0] - a[1];
		a[0] += a[1];
		for (j = 2; j < n; j += 2) {
			a[j - 1] = a[j] - a[j + 1];
			a[j] += a[j + 1];
		}
		a[n - 1] = xr;
	}
}


void ddst(int n, int isgn, double *a)
{
	void bitrv2(int n, double *a);
	void cftfsub(int n, double *a);
	void cftbsub(int n, double *a);
	void rftfsub(int n, double *a);
	void rftbsub(int n, double *a);
	void dstsub(int n, double *a);
	void dstsub4(int n, double *a);
	int j;
	double xr;

	if (isgn < 0) {
		xr = a[n - 1];
		for (j = n - 2; j >= 2; j -= 2) {
			a[j + 1] = -a[j] - a[j - 1];
			a[j] -= a[j - 1];
		}
		a[1] = a[0] + xr;
		a[0] -= xr;
		if (n > 4) {
			rftbsub(n, a);
			bitrv2(n, a);
			cftbsub(n, a);
		}
		else if (n == 4) {
			cftfsub(n, a);
		}
	}
	if (n > 4) {
		dstsub(n, a);
	}
	else {
		dstsub4(n, a);
	}
	if (isgn >= 0) {
		if (n > 4) {
			bitrv2(n, a);
			cftfsub(n, a);
			rftfsub(n, a);
		}
		else if (n == 4) {
			cftfsub(n, a);
		}
		xr = a[0] - a[1];
		a[0] += a[1];
		for (j = 2; j < n; j += 2) {
			a[j - 1] = -a[j] - a[j + 1];
			a[j] -= a[j + 1];
		}
		a[n - 1] = -xr;
	}
}


void dfct(int n, double *a)
{
	void ddct(int n, int isgn, double *a);
	void bitrv1(int n, double *a);
	int j, k, m, mh;
	double xr, xi, yr, yi, an;

	m = n >> 1;
	for (j = 0; j < m; j++) {
		k = n - j;
		xr = a[j] + a[k];
		a[j] -= a[k];
		a[k] = xr;
	}
	an = a[n];
	while (m >= 2) {
		ddct(m, 1, a);
		bitrv1(m, a);
		mh = m >> 1;
		xi = a[m];
		a[m] = a[0];
		a[0] = an - xi;
		an += xi;
		for (j = 1; j < mh; j++) {
			k = m - j;
			xr = a[m + k];
			xi = a[m + j];
			yr = a[j];
			yi = a[k];
			a[m + j] = yr;
			a[m + k] = yi;
			a[j] = xr - xi;
			a[k] = xr + xi;
		}
		xr = a[mh];
		a[mh] = a[m + mh];
		a[m + mh] = xr;
		m = mh;
	}
	xi = a[1];
	a[1] = a[0];
	a[0] = an + xi;
	a[n] = an - xi;
	bitrv1(n, a);
}


void dfst(int n, double *a)
{
	void ddst(int n, int isgn, double *a);
	void bitrv1(int n, double *a);
	int j, k, m, mh;
	double xr, xi, yr, yi;

	m = n >> 1;
	for (j = 1; j < m; j++) {
		k = n - j;
		xr = a[j] - a[k];
		a[j] += a[k];
		a[k] = xr;
	}
	a[0] = a[m];
	while (m >= 2) {
		ddst(m, 1, a);
		bitrv1(m, a);
		mh = m >> 1;
		for (j = 1; j < mh; j++) {
			k = m - j;
			xr = a[m + k];
			xi = a[m + j];
			yr = a[j];
			yi = a[k];
			a[m + j] = yr;
			a[m + k] = yi;
			a[j] = xr + xi;
			a[k] = xr - xi;
		}
		a[m] = a[0];
		a[0] = a[m + mh];
		a[m + mh] = a[mh];
		m = mh;
	}
	a[1] = a[0];
	a[0] = 0;
	bitrv1(n, a);
}


/* -------- child routines -------- */


#include <math.h>
#ifndef M_PI_2
#define M_PI_2      1.570796326794896619231321691639751442098584699687
#endif


#ifndef RDFT_LOOP_DIV  /* control of the RDFT's speed & tolerance */
#define RDFT_LOOP_DIV 64
#endif

#ifndef DCST_LOOP_DIV  /* control of the DCT,DST's speed & tolerance */
#define DCST_LOOP_DIV 64
#endif


void bitrv2(int n, double *a)
{
	int j0, k0, j1, k1, l, m, i, j, k;
	double xr, xi, yr, yi;

	l = n >> 2;
	m = 2;
	while (m < l) {
		l >>= 1;
		m <<= 1;
	}
	if (m == l) {
		j0 = 0;
		for (k0 = 0; k0 < m; k0 += 2) {
			k = k0;
			for (j = j0; j < j0 + k0; j += 2) {
				xr = a[j];
				xi = a[j + 1];
				yr = a[k];
				yi = a[k + 1];
				a[j] = yr;
				a[j + 1] = yi;
				a[k] = xr;
				a[k + 1] = xi;
				j1 = j + m;
				k1 = k + 2 * m;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += m;
				k1 -= m;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += m;
				k1 += 2 * m;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				for (i = n >> 1; i > (k ^= i); i >>= 1);
			}
			j1 = j0 + k0 + m;
			k1 = j1 + m;
			xr = a[j1];
			xi = a[j1 + 1];
			yr = a[k1];
			yi = a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
			for (i = n >> 1; i > (j0 ^= i); i >>= 1);
		}
	}
	else {
		j0 = 0;
		for (k0 = 2; k0 < m; k0 += 2) {
			for (i = n >> 1; i >(j0 ^= i); i >>= 1);
			k = k0;
			for (j = j0; j < j0 + k0; j += 2) {
				xr = a[j];
				xi = a[j + 1];
				yr = a[k];
				yi = a[k + 1];
				a[j] = yr;
				a[j + 1] = yi;
				a[k] = xr;
				a[k + 1] = xi;
				j1 = j + m;
				k1 = k + m;
				xr = a[j1];
				xi = a[j1 + 1];
				yr = a[k1];
				yi = a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				for (i = n >> 1; i > (k ^= i); i >>= 1);
			}
		}
	}
}


void bitrv2conj(int n, double *a)
{
	int j0, k0, j1, k1, l, m, i, j, k;
	double xr, xi, yr, yi;

	l = n >> 2;
	m = 2;
	while (m < l) {
		l >>= 1;
		m <<= 1;
	}
	if (m == l) {
		j0 = 0;
		for (k0 = 0; k0 < m; k0 += 2) {
			k = k0;
			for (j = j0; j < j0 + k0; j += 2) {
				xr = a[j];
				xi = -a[j + 1];
				yr = a[k];
				yi = -a[k + 1];
				a[j] = yr;
				a[j + 1] = yi;
				a[k] = xr;
				a[k + 1] = xi;
				j1 = j + m;
				k1 = k + 2 * m;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += m;
				k1 -= m;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				j1 += m;
				k1 += 2 * m;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				for (i = n >> 1; i > (k ^= i); i >>= 1);
			}
			k1 = j0 + k0;
			a[k1 + 1] = -a[k1 + 1];
			j1 = k1 + m;
			k1 = j1 + m;
			xr = a[j1];
			xi = -a[j1 + 1];
			yr = a[k1];
			yi = -a[k1 + 1];
			a[j1] = yr;
			a[j1 + 1] = yi;
			a[k1] = xr;
			a[k1 + 1] = xi;
			k1 += m;
			a[k1 + 1] = -a[k1 + 1];
			for (i = n >> 1; i > (j0 ^= i); i >>= 1);
		}
	}
	else {
		a[1] = -a[1];
		a[m + 1] = -a[m + 1];
		j0 = 0;
		for (k0 = 2; k0 < m; k0 += 2) {
			for (i = n >> 1; i >(j0 ^= i); i >>= 1);
			k = k0;
			for (j = j0; j < j0 + k0; j += 2) {
				xr = a[j];
				xi = -a[j + 1];
				yr = a[k];
				yi = -a[k + 1];
				a[j] = yr;
				a[j + 1] = yi;
				a[k] = xr;
				a[k + 1] = xi;
				j1 = j + m;
				k1 = k + m;
				xr = a[j1];
				xi = -a[j1 + 1];
				yr = a[k1];
				yi = -a[k1 + 1];
				a[j1] = yr;
				a[j1 + 1] = yi;
				a[k1] = xr;
				a[k1 + 1] = xi;
				for (i = n >> 1; i > (k ^= i); i >>= 1);
			}
			k1 = j0 + k0;
			a[k1 + 1] = -a[k1 + 1];
			a[k1 + m + 1] = -a[k1 + m + 1];
		}
	}
}


void bitrv1(int n, double *a)
{
	int j0, k0, j1, k1, l, m, i, j, k;
	double x;

	l = n >> 2;
	m = 1;
	while (m < l) {
		l >>= 1;
		m <<= 1;
	}
	if (m == l) {
		j0 = 0;
		for (k0 = 0; k0 < m; k0++) {
			k = k0;
			for (j = j0; j < j0 + k0; j++) {
				x = a[j];
				a[j] = a[k];
				a[k] = x;
				j1 = j + m;
				k1 = k + 2 * m;
				x = a[j1];
				a[j1] = a[k1];
				a[k1] = x;
				j1 += m;
				k1 -= m;
				x = a[j1];
				a[j1] = a[k1];
				a[k1] = x;
				j1 += m;
				k1 += 2 * m;
				x = a[j1];
				a[j1] = a[k1];
				a[k1] = x;
				for (i = n >> 1; i > (k ^= i); i >>= 1);
			}
			j1 = j0 + k0 + m;
			k1 = j1 + m;
			x = a[j1];
			a[j1] = a[k1];
			a[k1] = x;
			for (i = n >> 1; i > (j0 ^= i); i >>= 1);
		}
	}
	else {
		j0 = 0;
		for (k0 = 1; k0 < m; k0++) {
			for (i = n >> 1; i >(j0 ^= i); i >>= 1);
			k = k0;
			for (j = j0; j < j0 + k0; j++) {
				x = a[j];
				a[j] = a[k];
				a[k] = x;
				j1 = j + m;
				k1 = k + m;
				x = a[j1];
				a[j1] = a[k1];
				a[k1] = x;
				for (i = n >> 1; i >(k ^= i); i >>= 1);
			}
		}
	}
}


void cftfsub(int n, double *a)
{
	void cft1st(int n, double *a);
	void cftmdl(int n, int l, double *a);
	int j, j1, j2, j3, l;
	double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

	l = 2;
	if (n > 8) {
		cft1st(n, a);
		l = 8;
		while ((l << 2) < n) {
			cftmdl(n, l, a);
			l <<= 2;
		}
	}
	if ((l << 2) == n) {
		for (j = 0; j < l; j += 2) {
			j1 = j + l;
			j2 = j1 + l;
			j3 = j2 + l;
			x0r = a[j] + a[j1];
			x0i = a[j + 1] + a[j1 + 1];
			x1r = a[j] - a[j1];
			x1i = a[j + 1] - a[j1 + 1];
			x2r = a[j2] + a[j3];
			x2i = a[j2 + 1] + a[j3 + 1];
			x3r = a[j2] - a[j3];
			x3i = a[j2 + 1] - a[j3 + 1];
			a[j] = x0r + x2r;
			a[j + 1] = x0i + x2i;
			a[j2] = x0r - x2r;
			a[j2 + 1] = x0i - x2i;
			a[j1] = x1r - x3i;
			a[j1 + 1] = x1i + x3r;
			a[j3] = x1r + x3i;
			a[j3 + 1] = x1i - x3r;
		}
	}
	else {
		for (j = 0; j < l; j += 2) {
			j1 = j + l;
			x0r = a[j] - a[j1];
			x0i = a[j + 1] - a[j1 + 1];
			a[j] += a[j1];
			a[j + 1] += a[j1 + 1];
			a[j1] = x0r;
			a[j1 + 1] = x0i;
		}
	}
}


void cftbsub(int n, double *a)
{
	void cft1st(int n, double *a);
	void cftmdl(int n, int l, double *a);
	int j, j1, j2, j3, l;
	double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

	l = 2;
	if (n > 8) {
		cft1st(n, a);
		l = 8;
		while ((l << 2) < n) {
			cftmdl(n, l, a);
			l <<= 2;
		}
	}
	if ((l << 2) == n) {
		for (j = 0; j < l; j += 2) {
			j1 = j + l;
			j2 = j1 + l;
			j3 = j2 + l;
			x0r = a[j] + a[j1];
			x0i = -a[j + 1] - a[j1 + 1];
			x1r = a[j] - a[j1];
			x1i = -a[j + 1] + a[j1 + 1];
			x2r = a[j2] + a[j3];
			x2i = a[j2 + 1] + a[j3 + 1];
			x3r = a[j2] - a[j3];
			x3i = a[j2 + 1] - a[j3 + 1];
			a[j] = x0r + x2r;
			a[j + 1] = x0i - x2i;
			a[j2] = x0r - x2r;
			a[j2 + 1] = x0i + x2i;
			a[j1] = x1r - x3i;
			a[j1 + 1] = x1i - x3r;
			a[j3] = x1r + x3i;
			a[j3 + 1] = x1i + x3r;
		}
	}
	else {
		for (j = 0; j < l; j += 2) {
			j1 = j + l;
			x0r = a[j] - a[j1];
			x0i = -a[j + 1] + a[j1 + 1];
			a[j] += a[j1];
			a[j + 1] = -a[j + 1] - a[j1 + 1];
			a[j1] = x0r;
			a[j1 + 1] = x0i;
		}
	}
}


void cft1st(int n, double *a)
{
	int j, kj, kr;
	double ew, wn4r, wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
	double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

	x0r = a[0] + a[2];
	x0i = a[1] + a[3];
	x1r = a[0] - a[2];
	x1i = a[1] - a[3];
	x2r = a[4] + a[6];
	x2i = a[5] + a[7];
	x3r = a[4] - a[6];
	x3i = a[5] - a[7];
	a[0] = x0r + x2r;
	a[1] = x0i + x2i;
	a[4] = x0r - x2r;
	a[5] = x0i - x2i;
	a[2] = x1r - x3i;
	a[3] = x1i + x3r;
	a[6] = x1r + x3i;
	a[7] = x1i - x3r;
	wn4r = cos(M_PI_2 * 0.5);
	x0r = a[8] + a[10];
	x0i = a[9] + a[11];
	x1r = a[8] - a[10];
	x1i = a[9] - a[11];
	x2r = a[12] + a[14];
	x2i = a[13] + a[15];
	x3r = a[12] - a[14];
	x3i = a[13] - a[15];
	a[8] = x0r + x2r;
	a[9] = x0i + x2i;
	a[12] = x2i - x0i;
	a[13] = x0r - x2r;
	x0r = x1r - x3i;
	x0i = x1i + x3r;
	a[10] = wn4r * (x0r - x0i);
	a[11] = wn4r * (x0r + x0i);
	x0r = x3i + x1r;
	x0i = x3r - x1i;
	a[14] = wn4r * (x0i - x0r);
	a[15] = wn4r * (x0i + x0r);
	ew = M_PI_2 / n;
	kr = 0;
	for (j = 16; j < n; j += 16) {
		for (kj = n >> 2; kj >(kr ^= kj); kj >>= 1);
		wk1r = cos(ew * kr);
		wk1i = sin(ew * kr);
		wk2r = 1 - 2 * wk1i * wk1i;
		wk2i = 2 * wk1i * wk1r;
		wk3r = wk1r - 2 * wk2i * wk1i;
		wk3i = 2 * wk2i * wk1r - wk1i;
		x0r = a[j] + a[j + 2];
		x0i = a[j + 1] + a[j + 3];
		x1r = a[j] - a[j + 2];
		x1i = a[j + 1] - a[j + 3];
		x2r = a[j + 4] + a[j + 6];
		x2i = a[j + 5] + a[j + 7];
		x3r = a[j + 4] - a[j + 6];
		x3i = a[j + 5] - a[j + 7];
		a[j] = x0r + x2r;
		a[j + 1] = x0i + x2i;
		x0r -= x2r;
		x0i -= x2i;
		a[j + 4] = wk2r * x0r - wk2i * x0i;
		a[j + 5] = wk2r * x0i + wk2i * x0r;
		x0r = x1r - x3i;
		x0i = x1i + x3r;
		a[j + 2] = wk1r * x0r - wk1i * x0i;
		a[j + 3] = wk1r * x0i + wk1i * x0r;
		x0r = x1r + x3i;
		x0i = x1i - x3r;
		a[j + 6] = wk3r * x0r - wk3i * x0i;
		a[j + 7] = wk3r * x0i + wk3i * x0r;
		x0r = wn4r * (wk1r - wk1i);
		wk1i = wn4r * (wk1r + wk1i);
		wk1r = x0r;
		wk3r = wk1r - 2 * wk2r * wk1i;
		wk3i = 2 * wk2r * wk1r - wk1i;
		x0r = a[j + 8] + a[j + 10];
		x0i = a[j + 9] + a[j + 11];
		x1r = a[j + 8] - a[j + 10];
		x1i = a[j + 9] - a[j + 11];
		x2r = a[j + 12] + a[j + 14];
		x2i = a[j + 13] + a[j + 15];
		x3r = a[j + 12] - a[j + 14];
		x3i = a[j + 13] - a[j + 15];
		a[j + 8] = x0r + x2r;
		a[j + 9] = x0i + x2i;
		x0r -= x2r;
		x0i -= x2i;
		a[j + 12] = -wk2i * x0r - wk2r * x0i;
		a[j + 13] = -wk2i * x0i + wk2r * x0r;
		x0r = x1r - x3i;
		x0i = x1i + x3r;
		a[j + 10] = wk1r * x0r - wk1i * x0i;
		a[j + 11] = wk1r * x0i + wk1i * x0r;
		x0r = x1r + x3i;
		x0i = x1i - x3r;
		a[j + 14] = wk3r * x0r - wk3i * x0i;
		a[j + 15] = wk3r * x0i + wk3i * x0r;
	}
}


void cftmdl(int n, int l, double *a)
{
	int j, j1, j2, j3, k, kj, kr, m, m2;
	double ew, wn4r, wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
	double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

	m = l << 2;
	for (j = 0; j < l; j += 2) {
		j1 = j + l;
		j2 = j1 + l;
		j3 = j2 + l;
		x0r = a[j] + a[j1];
		x0i = a[j + 1] + a[j1 + 1];
		x1r = a[j] - a[j1];
		x1i = a[j + 1] - a[j1 + 1];
		x2r = a[j2] + a[j3];
		x2i = a[j2 + 1] + a[j3 + 1];
		x3r = a[j2] - a[j3];
		x3i = a[j2 + 1] - a[j3 + 1];
		a[j] = x0r + x2r;
		a[j + 1] = x0i + x2i;
		a[j2] = x0r - x2r;
		a[j2 + 1] = x0i - x2i;
		a[j1] = x1r - x3i;
		a[j1 + 1] = x1i + x3r;
		a[j3] = x1r + x3i;
		a[j3 + 1] = x1i - x3r;
	}
	wn4r = cos(M_PI_2 * 0.5);
	for (j = m; j < l + m; j += 2) {
		j1 = j + l;
		j2 = j1 + l;
		j3 = j2 + l;
		x0r = a[j] + a[j1];
		x0i = a[j + 1] + a[j1 + 1];
		x1r = a[j] - a[j1];
		x1i = a[j + 1] - a[j1 + 1];
		x2r = a[j2] + a[j3];
		x2i = a[j2 + 1] + a[j3 + 1];
		x3r = a[j2] - a[j3];
		x3i = a[j2 + 1] - a[j3 + 1];
		a[j] = x0r + x2r;
		a[j + 1] = x0i + x2i;
		a[j2] = x2i - x0i;
		a[j2 + 1] = x0r - x2r;
		x0r = x1r - x3i;
		x0i = x1i + x3r;
		a[j1] = wn4r * (x0r - x0i);
		a[j1 + 1] = wn4r * (x0r + x0i);
		x0r = x3i + x1r;
		x0i = x3r - x1i;
		a[j3] = wn4r * (x0i - x0r);
		a[j3 + 1] = wn4r * (x0i + x0r);
	}
	ew = M_PI_2 / n;
	kr = 0;
	m2 = 2 * m;
	for (k = m2; k < n; k += m2) {
		for (kj = n >> 2; kj >(kr ^= kj); kj >>= 1);
		wk1r = cos(ew * kr);
		wk1i = sin(ew * kr);
		wk2r = 1 - 2 * wk1i * wk1i;
		wk2i = 2 * wk1i * wk1r;
		wk3r = wk1r - 2 * wk2i * wk1i;
		wk3i = 2 * wk2i * wk1r - wk1i;
		for (j = k; j < l + k; j += 2) {
			j1 = j + l;
			j2 = j1 + l;
			j3 = j2 + l;
			x0r = a[j] + a[j1];
			x0i = a[j + 1] + a[j1 + 1];
			x1r = a[j] - a[j1];
			x1i = a[j + 1] - a[j1 + 1];
			x2r = a[j2] + a[j3];
			x2i = a[j2 + 1] + a[j3 + 1];
			x3r = a[j2] - a[j3];
			x3i = a[j2 + 1] - a[j3 + 1];
			a[j] = x0r + x2r;
			a[j + 1] = x0i + x2i;
			x0r -= x2r;
			x0i -= x2i;
			a[j2] = wk2r * x0r - wk2i * x0i;
			a[j2 + 1] = wk2r * x0i + wk2i * x0r;
			x0r = x1r - x3i;
			x0i = x1i + x3r;
			a[j1] = wk1r * x0r - wk1i * x0i;
			a[j1 + 1] = wk1r * x0i + wk1i * x0r;
			x0r = x1r + x3i;
			x0i = x1i - x3r;
			a[j3] = wk3r * x0r - wk3i * x0i;
			a[j3 + 1] = wk3r * x0i + wk3i * x0r;
		}
		x0r = wn4r * (wk1r - wk1i);
		wk1i = wn4r * (wk1r + wk1i);
		wk1r = x0r;
		wk3r = wk1r - 2 * wk2r * wk1i;
		wk3i = 2 * wk2r * wk1r - wk1i;
		for (j = k + m; j < l + (k + m); j += 2) {
			j1 = j + l;
			j2 = j1 + l;
			j3 = j2 + l;
			x0r = a[j] + a[j1];
			x0i = a[j + 1] + a[j1 + 1];
			x1r = a[j] - a[j1];
			x1i = a[j + 1] - a[j1 + 1];
			x2r = a[j2] + a[j3];
			x2i = a[j2 + 1] + a[j3 + 1];
			x3r = a[j2] - a[j3];
			x3i = a[j2 + 1] - a[j3 + 1];
			a[j] = x0r + x2r;
			a[j + 1] = x0i + x2i;
			x0r -= x2r;
			x0i -= x2i;
			a[j2] = -wk2i * x0r - wk2r * x0i;
			a[j2 + 1] = -wk2i * x0i + wk2r * x0r;
			x0r = x1r - x3i;
			x0i = x1i + x3r;
			a[j1] = wk1r * x0r - wk1i * x0i;
			a[j1 + 1] = wk1r * x0i + wk1i * x0r;
			x0r = x1r + x3i;
			x0i = x1i - x3r;
			a[j3] = wk3r * x0r - wk3i * x0i;
			a[j3 + 1] = wk3r * x0i + wk3i * x0r;
		}
	}
}


void rftfsub(int n, double *a)
{
	int i, j, k;
	double ec, w1r, w1i, wkr, wki, wdr, wdi, ss, xr, xi, yr, yi;

	ec = 2 * M_PI_2 / n;
	wkr = 0;
	wki = 0;
	wdi = cos(ec);
	wdr = sin(ec);
	wdi *= wdr;
	wdr *= wdr;
	w1r = 1 - 2 * wdr;
	w1i = 2 * wdi;
	ss = 2 * w1i;
	j = n >> 1;
	while (j > 4 * RDFT_LOOP_DIV) {
		for (i = 0; i < RDFT_LOOP_DIV; i++) {
			j -= 4;
			k = n - j;
			xr = a[j + 2] - a[k - 2];
			xi = a[j + 3] + a[k - 1];
			yr = wdr * xr - wdi * xi;
			yi = wdr * xi + wdi * xr;
			a[j + 2] -= yr;
			a[j + 3] -= yi;
			a[k - 2] += yr;
			a[k - 1] -= yi;
			wkr += ss * wdi;
			wki += ss * (0.5 - wdr);
			xr = a[j] - a[k];
			xi = a[j + 1] + a[k + 1];
			yr = wkr * xr - wki * xi;
			yi = wkr * xi + wki * xr;
			a[j] -= yr;
			a[j + 1] -= yi;
			a[k] += yr;
			a[k + 1] -= yi;
			wdr += ss * wki;
			wdi += ss * (0.5 - wkr);
		}
		wkr = 0.5 * sin(ec * j);
		wki = 0.5 * cos(ec * j);
		wdr = 0.5 - (wkr * w1r - wki * w1i);
		wdi = wkr * w1i + wki * w1r;
		wkr = 0.5 - wkr;
	}
	while (j > 4) {
		j -= 4;
		k = n - j;
		xr = a[j + 2] - a[k - 2];
		xi = a[j + 3] + a[k - 1];
		yr = wdr * xr - wdi * xi;
		yi = wdr * xi + wdi * xr;
		a[j + 2] -= yr;
		a[j + 3] -= yi;
		a[k - 2] += yr;
		a[k - 1] -= yi;
		wkr += ss * wdi;
		wki += ss * (0.5 - wdr);
		xr = a[j] - a[k];
		xi = a[j + 1] + a[k + 1];
		yr = wkr * xr - wki * xi;
		yi = wkr * xi + wki * xr;
		a[j] -= yr;
		a[j + 1] -= yi;
		a[k] += yr;
		a[k + 1] -= yi;
		wdr += ss * wki;
		wdi += ss * (0.5 - wkr);
	}
	xr = a[2] - a[n - 2];
	xi = a[3] + a[n - 1];
	yr = wdr * xr - wdi * xi;
	yi = wdr * xi + wdi * xr;
	a[2] -= yr;
	a[3] -= yi;
	a[n - 2] += yr;
	a[n - 1] -= yi;
}


void rftbsub(int n, double *a)
{
	int i, j, k;
	double ec, w1r, w1i, wkr, wki, wdr, wdi, ss, xr, xi, yr, yi;

	ec = 2 * M_PI_2 / n;
	wkr = 0;
	wki = 0;
	wdi = cos(ec);
	wdr = sin(ec);
	wdi *= wdr;
	wdr *= wdr;
	w1r = 1 - 2 * wdr;
	w1i = 2 * wdi;
	ss = 2 * w1i;
	j = n >> 1;
	a[j + 1] = -a[j + 1];
	while (j > 4 * RDFT_LOOP_DIV) {
		for (i = 0; i < RDFT_LOOP_DIV; i++) {
			j -= 4;
			k = n - j;
			xr = a[j + 2] - a[k - 2];
			xi = a[j + 3] + a[k - 1];
			yr = wdr * xr + wdi * xi;
			yi = wdr * xi - wdi * xr;
			a[j + 2] -= yr;
			a[j + 3] = yi - a[j + 3];
			a[k - 2] += yr;
			a[k - 1] = yi - a[k - 1];
			wkr += ss * wdi;
			wki += ss * (0.5 - wdr);
			xr = a[j] - a[k];
			xi = a[j + 1] + a[k + 1];
			yr = wkr * xr + wki * xi;
			yi = wkr * xi - wki * xr;
			a[j] -= yr;
			a[j + 1] = yi - a[j + 1];
			a[k] += yr;
			a[k + 1] = yi - a[k + 1];
			wdr += ss * wki;
			wdi += ss * (0.5 - wkr);
		}
		wkr = 0.5 * sin(ec * j);
		wki = 0.5 * cos(ec * j);
		wdr = 0.5 - (wkr * w1r - wki * w1i);
		wdi = wkr * w1i + wki * w1r;
		wkr = 0.5 - wkr;
	}
	while (j > 4) {
		j -= 4;
		k = n - j;
		xr = a[j + 2] - a[k - 2];
		xi = a[j + 3] + a[k - 1];
		yr = wdr * xr + wdi * xi;
		yi = wdr * xi - wdi * xr;
		a[j + 2] -= yr;
		a[j + 3] = yi - a[j + 3];
		a[k - 2] += yr;
		a[k - 1] = yi - a[k - 1];
		wkr += ss * wdi;
		wki += ss * (0.5 - wdr);
		xr = a[j] - a[k];
		xi = a[j + 1] + a[k + 1];
		yr = wkr * xr + wki * xi;
		yi = wkr * xi - wki * xr;
		a[j] -= yr;
		a[j + 1] = yi - a[j + 1];
		a[k] += yr;
		a[k + 1] = yi - a[k + 1];
		wdr += ss * wki;
		wdi += ss * (0.5 - wkr);
	}
	xr = a[2] - a[n - 2];
	xi = a[3] + a[n - 1];
	yr = wdr * xr + wdi * xi;
	yi = wdr * xi - wdi * xr;
	a[2] -= yr;
	a[3] = yi - a[3];
	a[n - 2] += yr;
	a[n - 1] = yi - a[n - 1];
	a[1] = -a[1];
}


void dctsub(int n, double *a)
{
	int i, j, k, m;
	double ec, w1r, w1i, wkr, wki, wdr, wdi, ss, xr, xi, yr, yi;

	ec = M_PI_2 / n;
	wkr = 0.5;
	wki = 0.5;
	w1r = cos(ec);
	w1i = sin(ec);
	wdr = 0.5 * (w1r - w1i);
	wdi = 0.5 * (w1r + w1i);
	ss = 2 * w1i;
	m = n >> 1;
	j = 0;
	while (j < m - 2 * DCST_LOOP_DIV) {
		for (i = 0; i < DCST_LOOP_DIV; i++) {
			j += 2;
			k = n - j;
			xr = wdi * a[j - 1] - wdr * a[k + 1];
			xi = wdr * a[j - 1] + wdi * a[k + 1];
			wkr -= ss * wdi;
			wki += ss * wdr;
			yr = wki * a[j] - wkr * a[k];
			yi = wkr * a[j] + wki * a[k];
			wdr -= ss * wki;
			wdi += ss * wkr;
			a[k + 1] = xr;
			a[k] = yr;
			a[j - 1] = xi;
			a[j] = yi;
		}
		wdr = cos(ec * j);
		wdi = sin(ec * j);
		wkr = 0.5 * (wdr - wdi);
		wki = 0.5 * (wdr + wdi);
		wdr = wkr * w1r - wki * w1i;
		wdi = wkr * w1i + wki * w1r;
	}
	while (j < m - 2) {
		j += 2;
		k = n - j;
		xr = wdi * a[j - 1] - wdr * a[k + 1];
		xi = wdr * a[j - 1] + wdi * a[k + 1];
		wkr -= ss * wdi;
		wki += ss * wdr;
		yr = wki * a[j] - wkr * a[k];
		yi = wkr * a[j] + wki * a[k];
		wdr -= ss * wki;
		wdi += ss * wkr;
		a[k + 1] = xr;
		a[k] = yr;
		a[j - 1] = xi;
		a[j] = yi;
	}
	xr = wdi * a[m - 1] - wdr * a[m + 1];
	a[m - 1] = wdr * a[m - 1] + wdi * a[m + 1];
	a[m + 1] = xr;
	a[m] *= wki + ss * wdr;
}


void dstsub(int n, double *a)
{
	int i, j, k, m;
	double ec, w1r, w1i, wkr, wki, wdr, wdi, ss, xr, xi, yr, yi;

	ec = M_PI_2 / n;
	wkr = 0.5;
	wki = 0.5;
	w1r = cos(ec);
	w1i = sin(ec);
	wdr = 0.5 * (w1r - w1i);
	wdi = 0.5 * (w1r + w1i);
	ss = 2 * w1i;
	m = n >> 1;
	j = 0;
	while (j < m - 2 * DCST_LOOP_DIV) {
		for (i = 0; i < DCST_LOOP_DIV; i++) {
			j += 2;
			k = n - j;
			xr = wdi * a[k + 1] - wdr * a[j - 1];
			xi = wdr * a[k + 1] + wdi * a[j - 1];
			wkr -= ss * wdi;
			wki += ss * wdr;
			yr = wki * a[k] - wkr * a[j];
			yi = wkr * a[k] + wki * a[j];
			wdr -= ss * wki;
			wdi += ss * wkr;
			a[j - 1] = xr;
			a[j] = yr;
			a[k + 1] = xi;
			a[k] = yi;
		}
		wdr = cos(ec * j);
		wdi = sin(ec * j);
		wkr = 0.5 * (wdr - wdi);
		wki = 0.5 * (wdr + wdi);
		wdr = wkr * w1r - wki * w1i;
		wdi = wkr * w1i + wki * w1r;
	}
	while (j < m - 2) {
		j += 2;
		k = n - j;
		xr = wdi * a[k + 1] - wdr * a[j - 1];
		xi = wdr * a[k + 1] + wdi * a[j - 1];
		wkr -= ss * wdi;
		wki += ss * wdr;
		yr = wki * a[k] - wkr * a[j];
		yi = wkr * a[k] + wki * a[j];
		wdr -= ss * wki;
		wdi += ss * wkr;
		a[j - 1] = xr;
		a[j] = yr;
		a[k + 1] = xi;
		a[k] = yi;
	}
	xr = wdi * a[m + 1] - wdr * a[m - 1];
	a[m + 1] = wdr * a[m + 1] + wdi * a[m - 1];
	a[m - 1] = xr;
	a[m] *= wki + ss * wdr;
}


void dctsub4(int n, double *a)
{
	int m;
	double wki, wdr, wdi, xr;

	wki = cos(M_PI_2 * 0.5);
	m = n >> 1;
	if (m == 2) {
		wdr = wki * sin(M_PI_2 * 0.25);
		wdi = wki * cos(M_PI_2 * 0.25);
		xr = wdi * a[1] - wdr * a[3];
		a[1] = wdr * a[1] + wdi * a[3];
		a[3] = xr;
	}
	a[m] *= wki;
}


void dstsub4(int n, double *a)
{
	int m;
	double wki, wdr, wdi, xr;

	wki = cos(M_PI_2 * 0.5);
	m = n >> 1;
	if (m == 2) {
		wdr = wki * sin(M_PI_2 * 0.25);
		wdi = wki * cos(M_PI_2 * 0.25);
		xr = wdi * a[3] - wdr * a[1];
		a[3] = wdr * a[3] + wdi * a[1];
		a[1] = xr;
	}
	a[m] *= wki;
}

//########################################################3
//########################################################3
//  Filter 처리 루틴
//########################################################3
//########################################################3


/************************************************************************

fir_filter_array - Filter float array of data

Requires FILTER structure for coefficients.
Input and output arrays are of equal length
and are allocated by the calller.

void fir_filter_array(float *in,float *out,int in_len,FILTER *fir)

in            float pointer to input array
out           float pointer to output array
in_len        length of input and output arrays
FILTER *fir   pointer to FILTER structure

*************************************************************************/

void fir_filter_array(double *in, double *out, int in_len, FILTER *fir)
{
	unsigned int i, j, coef_len2, acc_length;
	double acc;
	double *in_ptr, *data_ptr, *coef_start, *coef_ptr, *in_end;

	/* set up for coefficients */
	coef_start = fir->coef;
	coef_len2 = (fir->length + 1) / 2;

	/* set up input data pointers */
	in_end = in + in_len - 1;
	in_ptr = in + coef_len2 - 1;

	/* initial value of accumulation length for startup */
	acc_length = coef_len2;

	for (i = 0; i < in_len; i++) {

		/* set up pointers for accumulation */

		data_ptr = in_ptr;
		coef_ptr = coef_start;

		/* do accumulation and write result */

		acc = (*coef_ptr++) * (*data_ptr--);
		for (j = 1; j < acc_length; j++)
			acc += (*coef_ptr++) * (*data_ptr--);
		*out++ = acc;

		/* check for end case */

		if (in_ptr == in_end) {
			acc_length--;       /* one shorter each time */
			coef_start++;       /* next coefficient each time */
		}

		/* if not at end, then check for startup, add to input pointer */

		else {
			if (acc_length < fir->length) acc_length++;
			in_ptr++;
		}
	}
}

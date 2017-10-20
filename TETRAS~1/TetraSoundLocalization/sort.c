#include "sort.h"
int sort(int *x, int len)
{
	int i, j, imsi, curr_i = 0;
	int temp[125] = { 0, };

	for (i = 0; i<len; i++)
		for (j = i + 1; j<len; j++)
		{
			if (x[i] >= x[j])
			{
				imsi = x[i];
				x[i] = x[j];
				x[j] = imsi;
			}
		}

	temp[0] = x[0];
	curr_i++;
	for (i = 0; i<len; i++) {
		if (x[i] != x[i + 1]) {
			temp[curr_i] = x[i + 1];
			curr_i++;
		}
	}
	curr_i--;
	for (i = 0; i<len; i++) {
		if (i <= curr_i) {
			x[i] = temp[i];
		}
		else {
			x[i] = 0;
		}

	}
	return curr_i - 1;
}


double sortd(double *x, int len)
{
	int i, j, imsi, curr_i = 0;
	double temp[125] = { 0, };

	for (i = 0; i < len; i++)
		for (j = i + 1; j<len; j++)
		{
			if (x[i] >= x[j])
			{
				imsi = x[i];
				x[i] = x[j];
				x[j] = imsi;
			}
		}

	temp[0] = x[0];
	curr_i++;
	for (i = 0; i<len; i++) {
		if (x[i] != x[i + 1]) {
			temp[curr_i] = x[i + 1];
			curr_i++;
		}
	}
	curr_i--;
	for (i = 0; i<len; i++)
	{
		if (i <= curr_i)
			x[i] = temp[i];

		else
			x[i] = 0;
	}

	return curr_i - 1;
}
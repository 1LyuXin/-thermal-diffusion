#include <stdio.h>
#include <cmath>
#include <time.h>

int main(int argc, char* argv[])
{
	int N = 100000;
	float boundTemp = 10;
	float v = 0;
	float vv = 0.1;
	float* cur_Temp = new float [N];
	float* next_temp = new float [N];
	float* diff = new float [N];
	for (int i = 0; i < N; i++)
	{
		cur_Temp[i] = boundTemp;
	}
	for (int i = 1; i < N - 1; i++)
	{
		cur_Temp[i] = 0;
	}
	int iter = 0;
	clock_t t1 = clock();
	while (true)
	{
		v = 0;
		for (int i = 1; i < N - 1; i++)
		{
				next_temp[i] = (cur_Temp[i - 1] + cur_Temp[i + 1]) / 4;
				diff[i] = next_temp[i];
				v = (cur_Temp[i] - diff[i]) * (cur_Temp[i] - diff[i]) + v;
		}
		for (int i = 1; i < N - 1; i++)
		{
				cur_Temp[i] = diff[i];
		}
		iter++;
		if (v < vv) {
			break;
		}
	}
	clock_t t2 = clock();
	printf("N = %d, spend time %lf second\n", N, (double)(t2 - t1) / CLOCKS_PER_SEC);

	return 0;
}
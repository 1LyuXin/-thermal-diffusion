#include <stdio.h>
#include <cmath>
#include <time.h>

int main(int argc, char* argv[])
{
	int N = 1000;
	float boundTemp = 10;
	float v = 0;
	float vv = 0.1;
	float** cur_Temp = new float* [N];
	float** next_temp = new float* [N];
	float** diff = new float* [N];
	for (int i = 0; i < N; i++) {
		cur_Temp[i] = new float[N];
		next_temp[i] = new float[N];
		diff[i] = new float[N];
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cur_Temp[i][j] = boundTemp;
		}
	}
	for (int i = 1; i < N - 1; i++)
	{
		for (int j = 1; j < N - 1; j++)
		{
			cur_Temp[i][j] = 0;
		}
	}
	int iter = 0;
	clock_t t1 = clock();
	while (true)
	{
		v = 0;
		for (int i = 1; i < N - 1; i++)
		{
			for (int j = 1; j < N - 1; j++)
			{
				next_temp[i][j] = (cur_Temp[i - 1][j] + cur_Temp[i + 1][j] + cur_Temp[i][j - 1] + cur_Temp[i][j + 1]) / 4;
				diff[i][j] = next_temp[i][j];
				v = (cur_Temp[i][j] - diff[i][j]) * (cur_Temp[i][j] - diff[i][j]) + v;
			}
		}
		for (int i = 1; i < N - 1; i++)
		{
			for (int j = 1; j < N - 1; j++)
			{
				cur_Temp[i][j] = diff[i][j];
			}
		}
		iter++;
		if (v < vv){
			break;
		}
	}
	clock_t t2 = clock();
	printf("N = %d, spend time %lf second\n", N,(double)(t2 - t1) / CLOCKS_PER_SEC);

	return 0;
}
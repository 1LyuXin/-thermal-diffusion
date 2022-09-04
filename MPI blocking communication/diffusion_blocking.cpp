#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define DIFF 0.1
#define INIT_TEMP_BOUND 10
#define ITERS_MAX 10000
#define SIZE_MAX (500+2)


int block(float* cur_temp, int rank, int process_num, int N, int M, int iter) {
	int stop = 0;
	float* next_temp = (float*)malloc(sizeof(float) * N * M);
	MPI_Status status;
	int* pA;
	int* pB;
	int chunk = 1;
	float  stopTest = 0;
	if (rank > 0 && rank < process_num - 1) {
		MPI_Sendrecv(&cur_temp[N * (M - 2)], N, MPI_FLOAT, rank + 1, 0, &cur_temp[N * 0], N, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&cur_temp[N * 1], N, MPI_FLOAT, rank - 1, 1, &cur_temp[N * (M - 1)], N, MPI_FLOAT, rank + 1, 1, MPI_COMM_WORLD, &status);
	}
	else if (rank == 0 && rank < process_num - 1) {
		MPI_Sendrecv(&cur_temp[N * (M - 2)], N, MPI_FLOAT, rank + 1, 0, &cur_temp[N * (M - 1)], N, MPI_FLOAT, rank + 1, 1, MPI_COMM_WORLD, &status);
	}
	else if (rank > 0 && rank == process_num - 1) {
		MPI_Sendrecv(&cur_temp[N * 1], N, MPI_FLOAT, rank - 1, 1, &cur_temp[N * 0], N, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &status);
	}
	for (int i = 1; i < M - 1; i++) {
		for (int j = 1; j < N - 1; j++) {
			next_temp[N * i + j] = 0.25f * (cur_temp[N * (i - 1) + j] + cur_temp[N * (i + 1) + j] + cur_temp[N * i + j - 1] + cur_temp[N * i + j + 1]);
		}
	}
	pA = (int*)malloc(process_num * chunk * sizeof(int));
	if (!pA) {
		perror("can't allocate send buffer for stop test");
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}
	pB = (int*)malloc(process_num * chunk * sizeof(int));
	if (!pB) {
		perror("can't allocate receive buffer for stop test");
		free(pA);
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}
	for (int i = 0; i < process_num; i++) {
		for (int j = 0; j < chunk; j++) {
			pA[i * chunk + j] = 0;
			pB[i * chunk + j] = 0;
		}
	}

	if (iter % 100 == 0) {
		for (int i = 1; i < M - 1; i++) {
			for (int j = 1; j < N - 1; j++) {
				stopTest += (cur_temp[N * i + j] - next_temp[N * i + j]) * (cur_temp[N * i + j] - next_temp[N * i + j]);
			}
		}
		if (stopTest < DIFF) {
			for (int i = 0; i < process_num; i++) {
				for (int j = 0; j < chunk; j++) {
					pA[i * chunk + j] = 1;
				}
			}
		}
		stopTest = 0;
	}

	MPI_Alltoall(pA, chunk, MPI_INT, pB, chunk, MPI_INT, MPI_COMM_WORLD);

	for (int i = 0; i < process_num; i++) {
		for (int j = 0; j < chunk; j++) {
			if (pB[i * chunk + j] == 1) {
				stop = 1;
			}
		}
	}
	free(pA);
	free(pB);
	for (int i = 1; i < M - 1; i++) {
		for (int j = 1; j < N - 1; j++) {
			cur_temp[N * i + j] = next_temp[N * i + j];
		}
	}
	return stop;
}

int main(int argc, char* argv[]) {
	int M, N, rank, process_num;
	int step;
	float* cur_temp;
	double tStart, tEnd;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &process_num);


	N = SIZE_MAX;
	M = (N - 2) / process_num + 2;
	cur_temp = (float*)malloc(sizeof(float) * N * M);

	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			cur_temp[N * i + j] = 0;

	if (rank == 0) {
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {
				if ((i == 0) || (j == 0) || (j == N - 1)) {
					cur_temp[N * i + j] = INIT_TEMP_BOUND;
				}
			}
		}
	}
	if (rank > 0 && rank < process_num - 1) {
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {
				if ((j == 0) || (j == N - 1)) {
					cur_temp[N * i + j] = INIT_TEMP_BOUND;
				}
			}
		}
	}


	if (rank == process_num - 1) {
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {
				if ((i == M - 1) || (j == 0) || (j == N - 1)) {
					cur_temp[N * i + j] = INIT_TEMP_BOUND;
				}
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);


	tStart = MPI_Wtime();

	for (step = 1; step <= ITERS_MAX; step++) {
		if (block(cur_temp, rank, process_num, N, M, step)) {
			break;
		}
	}

	tEnd = MPI_Wtime();
	if (rank == 0) {
		printf("the size N: %d, spend time %lf seconds\n", N - 2, tEnd - tStart);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	free(cur_temp);

	MPI_Finalize();

	return 0;
}


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>


#define DIFF 0.1
#define INIT_TEMP_BOUND 10
#define ITERS_MAX 10000
#define SIZE_MAX (1000+2)

int no_block(float* cur_temp, int rank, int procss_num, int N, int M, int step) {
	int stop = 0;
	float* next_temp = (float*)malloc(sizeof(float) * N * M);
	int dntag = 2 * step;
	int uptag = 2 * step + 1;
	MPI_Request sndreq[2], rcvreq[2];
	MPI_Status stat[2];
	int* pA;
	int* pB;
	int chunk = 1;
	float  stopTest = 0;
	if (rank < procss_num - 1)
	{
		MPI_Isend(&cur_temp[N * (M - 2)], N, MPI_FLOAT, rank + 1, dntag, MPI_COMM_WORLD, &sndreq[0]);
		MPI_Irecv(&cur_temp[N * (M - 1)], N, MPI_FLOAT, rank + 1, uptag, MPI_COMM_WORLD, &rcvreq[0]);
	}

	if (rank > 0)
	{
		MPI_Isend(&cur_temp[N * 1], N, MPI_FLOAT, rank - 1, uptag, MPI_COMM_WORLD, &sndreq[1]);
		MPI_Irecv(&cur_temp[N * 0], N, MPI_FLOAT, rank - 1, dntag, MPI_COMM_WORLD, &rcvreq[1]);
	}
	for (int i = 2; i < M - 2; i++)
		for (int j = 1; j < N - 1; j++)
			next_temp[N * i + j] = 0.25*(cur_temp[N * (i - 1) + j] + cur_temp[N * (i + 1) + j] + cur_temp[N * i + j - 1] + cur_temp[N * i + j + 1]);
	if (procss_num > 1)
	{
		if (rank == 0)
			MPI_Wait(&rcvreq[0], stat);
		else if (rank == procss_num - 1)
			MPI_Wait(&rcvreq[1], stat);
		else
			MPI_Waitall(2, rcvreq, stat);
	}
	for (int j = 1; j < N - 1; j++)
	{
		next_temp[N * 1 + j] = 0.25f * (cur_temp[N * 0 + j] + cur_temp[N * 2 + j] + cur_temp[N * 1 + j - 1] + cur_temp[N * 1 + j + 1]);
		next_temp[N * (M - 2) + j] = 0.25f * (cur_temp[N * (M - 3) + j] + cur_temp[N * (M - 1) + j] + cur_temp[N * (M - 2) + j - 1] + cur_temp[N * (M - 2) + j + 1]);
	}
	if (procss_num > 1)
	{
		if (rank == 0)
			MPI_Wait(&sndreq[0], stat);
		else if (rank == procss_num - 1)
			MPI_Wait(&sndreq[1], stat);
		else
			MPI_Waitall(2, sndreq, stat);
	}
	pA = (int*)malloc(procss_num * chunk * sizeof(int));
	if (!pA) {
		perror("can't allocate send buffer for stop test");
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}
	pB = (int*)malloc(procss_num * chunk * sizeof(int));
	if (!pB) {
		perror("can't allocate receive buffer for stop test");
		free(pA);
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}
	for (int i = 0; i < procss_num; i++) {
		for (int j = 0; j < chunk; j++) {
			pA[i * chunk + j] = 0;
			pB[i * chunk + j] = 0;
		}
	}

	if (step % 100 == 0) {
		for (int i = 1; i < M - 1; i++) {
			for (int j = 1; j < N - 1; j++) {
				stopTest += (cur_temp[N * i + j] - next_temp[N * i + j]) * (cur_temp[N * i + j] - next_temp[N * i + j]);
			}
		}
		if (stopTest < DIFF) {
			for (int i = 0; i < procss_num; i++) {
				for (int j = 0; j < chunk; j++) {
					pA[i * chunk + j] = 1;
				}
			}
		}
		stopTest = 0;
	}

	MPI_Alltoall(pA, chunk, MPI_INT, pB, chunk, MPI_INT, MPI_COMM_WORLD);

	for (int i = 0; i < procss_num; i++) {
		for (int j = 0; j < chunk; j++) {
			if (pB[i * chunk + j] == 1) {
				stop = 1;
			}
		}
	}
	free(pA);
	free(pB);
	for (int i = 1; i < M - 1; i++)
		for (int j = 1; j < N - 1; j++)
			cur_temp[N * i + j] = next_temp[N * i + j];

	return stop;
}

int main(int argc, char* argv[]) {
	int M, N, rank, process_num;
	int step;
	float* cur_temp;
	double t1, t2;
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
		for (int i = 0; i < M; i++)
			for (int j = 0; j < N; j++)
				if ((i == 0) || (j == 0) || (j == N - 1)) {
					cur_temp[N * i + j] = INIT_TEMP_BOUND;
				}
	}
	if (rank > 0 && rank < process_num - 1) {
		for (int i = 0; i < M; i++)
			for (int j = 0; j < N; j++)
				if ((j == 0) || (j == N - 1))
					cur_temp[N * i + j] = INIT_TEMP_BOUND;
	}
	if (rank == process_num - 1) {
		for (int i = 0; i < M; i++)
			for (int j = 0; j < N; j++)
				if ((i == M - 1) || (j == 0) || (j == N - 1))
					cur_temp[N * i + j] = INIT_TEMP_BOUND;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();
	for (step = 1; step <= ITERS_MAX; step++)
		if (no_block(cur_temp, rank, process_num, N, M, step))
			break;
	t2 = MPI_Wtime();

	if (rank == 0) {
		printf("the size N: %d, spend time %lf seconds\n", N - 2, t2 - t1);
	}


	MPI_Barrier(MPI_COMM_WORLD);

	free(cur_temp);

	MPI_Finalize();

	return 0;
}

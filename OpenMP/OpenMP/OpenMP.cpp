#include <iostream>
#include <omp.h>
#include <time.h>

#define N 2000
#define M 750
#define L 1500
using namespace std;

void matrix_mul(double** X, double** Y, double** Z) {
#pragma omp parallel for shared(X, Y, Z)
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < L; j++) {
			Z[i][j] = 0.0;
			for (int k = 0; k < M; k++)
				Z[i][j] += X[i][k] * Y[k][j];
		}
	}
}


int main()
{
	int i, j;

	double** x, ** y, ** z;
	double a1, a2;

	x = new double* [N];
	for (i = 0; i < N; i++)
		x[i] = new double[N];

	for (i = 0; i < N; i++)
		for (j = 0; j < M; j++)
			x[i][j] = i + j;

	y = new double* [L];
	for (i = 0; i < M; i++)
		y[i] = new double[L];

	for (i = 0; i < M; i++)
		for (j = 0; j < L; j++)
			y[i][j] = i + j;

	for (int l = 1; l <= omp_get_num_procs(); l++) {
		z = new double* [N];
		for (i = 0; i < N; i++)
			z[i] = new double[N];
		string threads = "threads";
		if (l == 1)
			threads = "thread";
		cout << l << " " << threads << endl;
		omp_set_num_threads(l);
		a1 = omp_get_wtime();
		matrix_mul(x, y, z);
		a2 = omp_get_wtime();
		cout << "Multiplication time " << a2 - a1 << " seconds" << endl;
		cout << "" << endl;
	}
}
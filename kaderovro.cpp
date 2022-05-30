﻿#include "kaderovro.h"

/**
 * Введение в дисциплину
 */
void kaderovro::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void kaderovro::lab2()
{
	for (int i = 0; i < N; i++) {

		int maxNum = i;

		for (int j = i + 1; j < N; j++) {
			if (abs(A[j][i]) > abs(A[maxNum][i])) {
				maxNum = j;
			}
		}

		if (maxNum != i) {
			for (int j = 0; j < N; j++) {
				swap(A[maxNum][j], A[i][j]);
				swap(b[maxNum], b[i]);
			}
		}

		for (int j = i + 1; j < N; j++) {
			double s = A[j][i] / A[i][i];
			for (int c = i; c < N; c++) {
				A[j][c] -= s * A[i][c];
			}
			b[j] -= (s * b[i]);
		}
	}
	for (int i = N - 1; i >= 0; i--) {
		for (int j = i + 1; j < N; j++) {
			b[i] -= A[i][j] * x[j];
		}
		x[i] = b[i] / A[i][i];
	}
		
	
}



/**
 * Метод прогонки
 */
void kaderovro::lab3()
{
	double* alpha = new double[N];
	double* beta = new double[N];
	double y;
	
	y = A[0][0];
	alpha[0] = -A[0][1] / y;
	beta[0] = b[0] / y;
	
	for (int i = 1; i < N; i++) {
		y = A[i][i] + A[i][i - 1] * alpha[i - 1];
		alpha[i] = -A[i][i+1] / y;
		beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / y;
	}
	x[N - 1] = beta[N - 1];
	
	for (int i = N - 2; i >= 0; i--) {
		x[i] = alpha[i] * x[i + 1] + beta[i];
	}
}



/**
 * Метод Холецкого
 */
void kaderovro::lab4()
{
	double** L = new double* [N];
	for (int i = 0; i < N; i++) {
		L[i] = new double[N];
		for (int j = 0; j < N; j++) L[i][j] = 0;
	}

	for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
			if (j == i) {
				for (int k = 0; k < i; k++) A[i][i] -= L[k][i] * L[k][i];
				L[i][i] = sqrt(A[i][i]);
				
			}
			else {
				for (int k = 0; k < j; k++) A[i][j] -= L[k][i] * L[k][j];
				L[i][j] = A[i][j] / L[i][i];
			}
		}
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < i; j++) b[i] -= L[j][i] * b[j];
		b[i] /= L[i][i];
	}

	for (int i = N - 1; i >= 0; i--) {
		for (int k = i + 1; k < N; k++) b[i] -= L[i][k] * x[k];
		x[i] = b[i] / L[i][i];
	}

	for (int i = 0; i < N; i++) delete[] L[i];
	delete[] L;
}



/**
 * Метод Якоби или Зейделя
 */
void kaderovro::lab5()
{
	double norm, sum, approx;
	double eps = 1.e-15;
	for (int i = 0; i < N; ++i) x[i] = b[i];

	do {
		norm = 0;
		for (int i = 0; i < N; i++) {
			sum = 0;
			for (int j = 0; j < N; j++) if (i != j) sum += A[i][j] * x[j];
			approx = (b[i] - sum) / A[i][i];
			if (norm < (abs(approx - x[i]))) norm = (abs(approx - x[i]));
			x[i] = approx;
		}
	} while (norm >= eps);
}



double multOnSkaclar(double* A, double* B, int N)
{
	double res = 0;
	for (int i = 0; i < N; ++i)
	{
		res += A[i] * B[i];
	}
	return res;
}


double* multOnVector(double** A, double* V, int N)
{
	double* res = new double[N];
	for (int i = 0; i < N; ++i)
	{
		res[i] = 0;
		for (int j = 0; j < N; ++j) res[i] += A[i][j] * V[j];
	}
	return res;
}

/**
 * Метод минимальных невязок
 */
void kaderovro::lab6()
{
	double norm = 1e-1;
	
	double eps = 1.e-15;
	double* L = new double [N];
	for (int i = 0; i < N; i++) L[i] = 0;

	for (int i = 0; i < N; ++i)
	{
		x[i] = 1e-2;
	}

	do
	{
		double* Ax = multOnVector(A, x, N);

		for (int i = 0; i < N; ++i)
		{
			L[i] = Ax[i] - b[i];
		}

		double* Ar = multOnVector(A, L, N);
		double tau = 0;

		tau = (multOnSkaclar(Ar, L, N) / multOnSkaclar(Ar, Ar, N));

		for (int i = 0; i < N; i++)
		{
			double prevX = x[i];
			x[i] = prevX - tau * L[i];

			if (fabs(x[i] - prevX) < norm)
			{
				norm = fabs(x[i] - prevX);
			}
		}

		delete[] Ar;
		delete[] Ax;

	} while (sqrt(norm) >= eps);

}



/**
 * Метод сопряженных градиентов
 */
void kaderovro::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void kaderovro::lab8()
{
	double eps = 1e-20;
	double** B = new double* [N];
	for (int i = 0; i < N; i++) {
		B[i] = new double[N];
	}

	while (true) {
		double norm = 0;
		int imax = 0;
		int jmax = 1;
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				if (abs(A[i][j]) > abs(A[imax][jmax])) {
					imax = i;
					jmax = j;
				}
				norm += A[i][j] * A[i][j];
			}
		}

		if (sqrt(norm) < eps) {
			break;
		}

		double fi = 0.5 * atan(2 * A[imax][jmax] / (A[imax][imax] - A[jmax][jmax]));

		for (int i = 0; i < N; i++) {
			B[i][imax] = A[i][imax] * cos(fi) + A[i][jmax] * sin(fi);
			B[i][jmax] = A[i][jmax] * cos(fi) - A[i][imax] * sin(fi);
		}

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (j != imax && j != jmax) {
					B[i][j] = A[i][j];
				}
			}
		}

		for (int j = 0; j < N; j++) {
			A[imax][j] = B[imax][j] * cos(fi) + B[jmax][j] * sin(fi);
			A[jmax][j] = B[jmax][j] * cos(fi) - B[imax][j] * sin(fi);
		}

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (i != imax && i != jmax) {
					A[i][j] = B[i][j];
				}
			}
		}
	}

	for (int i = 0; i < N; i++) {
		x[i] = A[i][i];
	}

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void kaderovro::lab9()
{
	double AbsMaxEigenvalue;
	double Eigenvalue = 0.;
	double eps = 1e-8;
	double* xPrev = new double[N];
	double* xNew = new double[N];

	for (int i = 0; i < N; ++i)
	{
		xPrev[i] = b[i];
	}

	do
	{
		AbsMaxEigenvalue = Eigenvalue;

		xNew = multOnVector(A, xPrev, N);

		double s1 = 0, s2 = 0;
		for (int i = 0; i < N; ++i)
		{
			s1 += xNew[i];
			s2 += xPrev[i];
		}

		Eigenvalue = s1 / s2;

		for (int i = 0; i < N; ++i)
		{
			xPrev[i] = xNew[i] / multOnSkaclar(xNew, xNew, N);
		}

	} while (fabs(AbsMaxEigenvalue - Eigenvalue) >= eps);

	delete[] xNew;
	delete[] xPrev;

	cout << "Absolute max eigenvalue = " << AbsMaxEigenvalue << endl;
}


std::string kaderovro::get_name()
{
  return "R.O. Kaderov";
}

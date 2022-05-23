#include "prokopenkods.h"

/**
 * Введение в дисциплину
 */
void prokopenkods::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void prokopenkods::lab2()
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
void prokopenkods::lab3()
{
	double* alpha = new double[N];
	double* beta = new double[N];
	double y;

	y = A[0][0];
	alpha[0] = -A[0][1] / y;
	beta[0] = b[0] / y;

	for (int i = 1; i < N; i++) {
		y = A[i][i] + A[i][i - 1] * alpha[i - 1];
		alpha[i] = -A[i][i + 1] / y;
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
void prokopenkods::lab4()
{
	double S[N][N];
	int D[N][N];
	for (int i = 0; i < N; i++) {

		double subSum = 0;
		for (int l = 0; l <= i - 1; l++) {
			subSum += S[l][i] * S[l][i] * D[l][l];
		}

		int sign = (A[i][i] - subSum) < 0;
		D[i][i] = (int)pow(-1, sign);

		S[i][i] = sqrt(abs(A[i][i] - subSum));
		for (int j = 0; j < i; j++) S[i][j] = 0.0;
		for (int j = i + 1; j < N; j++) {
			double subSum = 0;
			for (int l = 0; l <= i - 1; l++) subSum += S[l][j] * S[l][i] * D[l][l];
			S[i][j] = (A[i][j] - subSum) / S[i][i] * D[i][i];
		}

		//Обратный ход

		double y[N];
		y[0] = b[0] / S[0][0];
		for (int i = 1; i < N; i++) {
			double subSum = 0;
			for (int j = 0; j <= i - 1; j++) subSum += S[j][i] * y[j];
			y[i] = b[i] - subSum;
			y[i] /= S[i][i];
		}
		x[N - 1] = y[N - 1] / S[N - 1][N - 1];
		for (int i = N - 2; i >= 0; i--) {
			double subSum = 0;
			for (int k = i + 1; k <= N - 1; k++) subSum += S[i][k] * x[k];
			x[i] = y[i] - subSum;
			x[i] /= S[i][i];
		}

	}

}



/**
 * Метод Якоби или Зейделя
 */
void prokopenkods::lab5()
{
	double* xk = new double[N];
	double eps = 1e-20;
	double norm;
	for (int i = 0; i < N; i++) x[i] = b[i] / A[i][i];

	do {
		for (int i = 0; i < N; i++) {

			for (int j = 0; j < N; j++) xk[j] = x[j];

			double lowerSum = 0, upperSum = 0;

			for (int j = 0; j < i; j++) lowerSum += A[i][j] * xk[j];
			for (int j = i + 1; j < N; j++) upperSum += A[i][j] * xk[j];

			x[i] = 1 / A[i][i] * (b[i] - lowerSum - upperSum);

			if (i == 0) norm = abs(x[i] - xk[i]);
			if (abs(x[i] - xk[i]) > norm) norm = abs(x[i] - xk[i]);

		}

	} while (norm >= eps);
}



/**
 * Метод минимальных невязок
 */
void prokopenkods::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void prokopenkods::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void prokopenkods::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void prokopenkods::lab9()
{

}


std::string prokopenkods::get_name()
{
  return "D.S.Prokopenko";
}

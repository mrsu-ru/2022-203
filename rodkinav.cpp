#include "rodkinav.h"
/**
 * Введение в дисциплину
 */
 // Вспомогательный метод произведения матриц на вектор
double* Composition(double** A, double* x, int N) {
	double* result = new double[N];
	for (int i = 0; i < N; i++) {
		result[i] = 0;
		for (int j = 0; j < N; j++) {
			result[i] += A[i][j] * x[j];
		}
	}
	return result;
}
//Вспомогательный метод вычисления скалярного произведения двух векторов
double sca(double* x, double* y, int N) {
	double scal = 0;
	for (int i = 0; i < N; i++) scal += x[i] * y[i];
	return scal;
}
void rodkinav::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void rodkinav::lab2()
{
	int max;
	double sum = 0;
	//прямой ход
	for (int i = 0; i < N; i++)
	{
		max = i;
		double* vspom; double vspom_1;

		for (int j = i + 1; j < N; j++)
		{
			if (fabs(A[j][i]) > fabs(A[max][i]))
			{
				max = j;
			}
		}
		if (max != i)
		{
			vspom = A[i];
			A[i] = A[max];
			A[max] = vspom;

			vspom_1 = b[i];
			b[i] = b[max];
			b[max] = vspom_1;
		}

		for (int j = N - 1; j > i; j--)
		{
			A[i][j] /= A[i][i];
		}
		b[i] /= A[i][i];
		A[i][i] = 1;

		for (int j = i + 1; j < N; j++)
		{
			for (int k = N - 1; k > i; k--)
			{
				A[j][k] = A[j][k] - A[j][i] * A[i][k];
			}
			b[j] = b[j] - A[j][i] * b[i];
			A[j][i] = 0;
		}
	}

	//обратный ход
	x[N - 1] = b[N - 1];
	for (int i = N - 2; i > -1; i--)
	{
		for (int j = i + 1; j < N; j++)
		{
			sum += x[j] * A[i][j];
		}
		x[i] = b[i] - sum;
		sum = 0;
	}
}



/**
 * Метод прогонки
 */
void rodkinav::lab3()
{
	double* alfa = new double[N];
	double* beta = new double[N];
	alfa[0] = A[0][1] / -A[0][0];
	beta[0] = b[0] / A[0][0];
	for (int i = 1; i < N; i++) {
		alfa[i] = A[i][i + 1] / (-A[i][i] - A[i][i - 1] * alfa[i - 1]);
		beta[i] = (A[i][i - 1] * beta[i - 1] - b[i]) / (-A[i][i] - A[i][i - 1] * alfa[i - 1]);
	}
	x[N - 1] = beta[N - 1];
	for (int i = N - 2; i >= 0; i--) {
		x[i] = alfa[i] * x[i + 1] + beta[i];
	}
}



/**
 * Метод Холецкого
 */
void rodkinav::lab4()
{
	double S[N][N];
	int D[N][N];
	for (int i = 0; i < N; i++) 
	{

		double subSum = 0;
		for (int l = 0; l <= i - 1; l++) 
		{
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
		for (int i = 1; i < N; i++) 
		{
			double subSum = 0;
			for (int j = 0; j <= i - 1; j++) subSum += S[j][i] * y[j];
			y[i] = b[i] - subSum;
			y[i] /= S[i][i];
		}
		x[N - 1] = y[N - 1] / S[N - 1][N - 1];
		for (int i = N - 2; i >= 0; i--) 
		{
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
void rodkinav::lab5()
{
	double* xk = new double[N];
	double eps = 1e-20;
	double norma;
	//Начальное приближение
	for (int i = 0; i < N; i++) x[i] = b[i] / A[i][i];

	do {
		for (int i = 0; i < N; i++) 
		{
			for (int j = 0; j < N; j++) xk[j] = x[j];
			double lowerSum = 0, upperSum = 0;
			for (int j = 0; j < i; j++) lowerSum += A[i][j] * xk[j];
			for (int j = i + 1; j < N; j++) upperSum += A[i][j] * xk[j];
			x[i] = 1 / A[i][i] * (b[i] - lowerSum - upperSum);
			if (i == 0) norma = abs(x[i] - xk[i]);
			if (abs(x[i] - xk[i]) > norma) norma = abs(x[i] - xk[i]);

		}

	} while (norma >= eps);
}



/**
 * Метод минимальных невязок
 */
void rodkinav::lab6()
{
	double eps = 1e-18;
	double* rk = new double[N];
	double T = 0;
	for (int i = 0; i < N; i++) x[i] = 0.1;
	double norm;
	do {
		double* Axk = Composition(A, x, N);
		for (int i = 0; i < N; i++) rk[i] = Axk[i] - b[i];
		double* Ark = Composition(A, rk, N);
		T = sca(Ark, rk, N) / sca(Ark, Ark, N);
		norm = 0;
		for (int i = 0; i < N; i++) 
		{
			double check = x[i];
			x[i] = x[i] - T * rk[i];
			norm += (x[i] - check) * (x[i] - check);
		}
	} while (sqrt(norm) > eps);
}



/**
 * Метод сопряженных градиентов
 */
void rodkinav::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void rodkinav::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void rodkinav::lab9()
{

}


std::string rodkinav::get_name()
{
  return "A.V. Rodkin";
}

#include "kaderovro.h"

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



/**
 * Метод минимальных невязок
 */
void kaderovro::lab6()
{

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

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void kaderovro::lab9()
{

}


std::string kaderovro::get_name()
{
  return "R.O. Kaderov";
}

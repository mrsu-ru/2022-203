#include "kochetkovpa.h"

/**
 * Введение в дисциплину
 */
void kochetkovpa::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void kochetkovpa::lab2()
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
void kochetkovpa::lab3()
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
void kochetkovpa::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void kochetkovpa::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void kochetkovpa::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void kochetkovpa::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void kochetkovpa::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void kochetkovpa::lab9()
{

}


std::string kochetkovpa::get_name()
{
  return "P.A. Kochetkov";
}

#include "fedinda.h"
#include <iostream>
const double eps = 1.e-15;
/**
 * Введение в дисциплину
 */
using namespace std;

void fedinda::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void fedinda::lab2()
{
  for (int i = 0; i < N; i++) {
		int max = i;
		for (int j = i; j < N; j++) {
			if (A[j][i] > A[max][i]) {
				max = j;
			}
		}
		if (max != i) {
			swap(A[i], A[max]);
			swap(b[i], b[max]);
		}
		double mainElem = A[i][i];
		for (int j = i; j < N; j++) {
			A[i][j] /= mainElem;
		}
		b[i] /= mainElem;

		for (int j = 0; j < N; j++) {
			if (j != i) {
				mainElem = A[j][i];
				for (int k = 0; k < N;k++) {
					A[j][k] -= mainElem * A[i][k];
				}
				b[j] -= mainElem * b[i];
			}
		}
	}
  
	for (int i = 0; i < N;i++) {
		x[i] = b[i];
	}
}



/**
 * Метод прогонки
 */
void fedinda::lab3()
{
	double* alpha = new double[N];
	double* beta = new double[N];
	int i;
	alpha[0] = A[0][1] / A[0][0];
	beta[0] = b[0] / A[0][0];
	for (i = 1; i < N; i++) {
		alpha[i] = A[i][i + 1] / (A[i][i] - A[i][i - 1] * alpha[i - 1]);
		beta[i] = (b[i] - beta[i - 1] * A[1][0]) / (A[i][i] - A[i][i - 1] * alpha[i - 1]);
	}

	x[N - 1] = beta[N - 1];
	for (i = N - 2; i >= 0; i--) {
		x[i] = beta[i] - alpha[i] * x[i + 1];
	}
}



/**
 * Метод Холецкого
 */
void fedinda::lab4()
{
	const double lambda = 0.01;
	double prev[N];
	for (int i = 0; i < N; i++){
		prev[i] = 0;
	}

	double difference;
	do{
		for (int i = 0; i < N; i++){
			double sum = 0;
			for (int j = 0; j < N; j++){
				sum += A[i][j] * prev[j];
			}
			x[i] = prev[i] - lambda * (sum - b[i]);
		}
		difference = 0;
		for (int i = 0; i < N; i++){
			if (fabs(x[i] - prev[i]) > difference){
				difference = fabs(x[i] - prev[i]);
			}
		}
		for (int i = 0; i < N; i++){
			prev[i] = x[i];
		}
	} while (difference > eps);
}



/**
 * Метод Якоби или Зейделя
 */
void fedinda::lab5()
{
	double prev[N];
	for (int i = 0; i < N; i++)
	{
		prev[i] = 0;
	}

	double difference;
	do
	{
		for (int i = 0; i < N; i++)
		{
			double sum = 0;
			for (int j = 0; j < i; j++)
			{
				sum += A[i][j] * prev[j];
			}

			for (int j = i + 1; j < N; j++)
			{
				sum += A[i][j] * prev[j];
			}

			x[i] = (b[i] - sum) / A[i][i];
		}

		difference = 0;
		for (int i = 0; i < N; i++)
		{
			if (fabs(x[i] - prev[i]) > difference)
			{
				difference = fabs(x[i] - prev[i]);
			}
		}

		for (int i = 0; i < N; i++)
		{
			prev[i] = x[i];
		}

	} while (difference > eps);
}



/**
 * Метод минимальных невязок
 */
void fedinda::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void fedinda::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void fedinda::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void fedinda::lab9()
{

}


std::string fedinda::get_name()
{
  return "D.A. Fedin";
}

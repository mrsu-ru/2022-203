#include "timovkinayu.h"

/**
 * Введение в дисциплину
 */
void timovkinayu::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void timovkinayu::lab2()
{
    double p;
    int maxn;

    for (int k = 0; k < N - 1; k++) {
        maxn = k;
        for (int i = k + 1; i < N; i++)
            if (abs(A[i][k]) > abs(A[maxn][k])) maxn = i; 
        swap(A[maxn], A[k]);
        swap(b[maxn], b[k]);

        for (int i = k + 1; i < N; i++) {
            p = A[i][k] / A[k][k];
            for (int j = k; j < N; j++)
                A[i][j] -= p * A[k][j];
            b[i] -= p * b[k];
        }
    }

    for (int i = 0; i < N; i++) x[i] = b[i];
    for (int i = N - 1; i >= 0; i--) {
        for (int j = i + 1; j < N; j++)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
    for (int i = 0; i < N; i++) cout << x[i] << " ";
    cout << endl;
}



/**
 * Метод прогонки
 */
void timovkinayu::lab3()
{
	double* alpha = new double[N]; 
	double* betta = new double[N]; 

	alpha[0] = -A[0][1] / A[0][0];
	betta[0] = b[0] / A[0][0];

	for (int i = 1; i < N; i++) { 

		alpha[i] = A[i][i + 1] / (-A[i][i] - A[i][i - 1] * alpha[i - 1]);
		betta[i] = (-b[i] + A[i][i - 1] * betta[i - 1]) / (-A[i][i] - A[i][i - 1] * alpha[i - 1]);
	}


	x[N - 1] = betta[N - 1];
	for (int i = N - 2; i >= 0; i--) 
		x[i] = alpha[i] * x[i + 1] + betta[i];


	for (int i = 0; i < N; i++) cout << x[i] << " ";
	cout << endl;
	delete[] alpha;
	delete[] betta;
}



/**
 * Метод Холецкого
 */
void timovkinayu::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void timovkinayu::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void timovkinayu::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void timovkinayu::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void timovkinayu::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void timovkinayu::lab9()
{

}


std::string timovkinayu::get_name()
{
  return "A.YU. Timovkin";
}

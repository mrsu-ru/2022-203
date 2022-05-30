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

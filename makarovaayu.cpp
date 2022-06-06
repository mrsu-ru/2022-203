#include "makarovaayu.h"

/**
 * Введение в дисциплину
 */
void makarovaayu::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void makarovaayu::lab2()
{
    for (int i = 0; i < N; i++) {
        int maxLine = i;
        for (int j = i + 1; j < N; j++) {
            if (abs(A[j][i]) > abs(A[maxLine][i])) {
                maxLine = j;
            }
        }
        if (i != maxLine) {
            swap(A[maxLine], A[i]);
            swap(b[maxLine], b[i]);
        }
        for (int j = i + 1; j < N; j++) {
            double c = A[j][i] / A[i][i];
            for (int k = i; k < N; k++) {
                A[j][k] -= (c * A[i][k]);
            }
            b[j] -= (c * b[i]);
        }
    }

    for (int i = N - 1; i >= 0; i--) {
        for (int j = i + 1; j < N; j++) {
            b[i] -= (A[i][j] * x[j]);
        }
        x[i] = b[i] / A[i][i];
    }

}



/**
 * Метод прогонки
 */
void makarovaayu::lab3()
{
        double alpha[N - 1], beta[N];
        alpha[0] = - A[0][1] / A[0][0];
        beta[0] = b[0] / A[0][0];

        for (int i = 1; i < N; i++)
        {
            double y = A[i][i] + A[i][i - 1] * alpha[i - 1];
            alpha[i] = - A[i][i + 1] / y;
            beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / y;
        }

        x[N - 1] = beta[N - 1];
        for (int i = N - 2; i >= 0; i--)
        {
            x[i] = alpha[i] * x[i + 1] + beta[i];
        }
    }





/**
 * Метод Холецкого
 */
void makarovaayu::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void makarovaayu::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void makarovaayu::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void makarovaayu::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void makarovaayu::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void makarovaayu::lab9()
{

}


std::string makarovaayu::get_name()
{
  return "A.U.Makarova";
}

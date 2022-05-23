#include "denisovrv.h"

/**
 * Введение в дисциплину
 */
void denisovrv::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void denisovrv::lab2()
{
    for (int i = 0; i < N; i++) {
        int max = 0;
        for (int j = 0; j < N; j++) if (abs(A[j][i]) > abs(A[max][i])) max = j;
        swap(A[i], A[max]);
        swap(b[i], b[max]);
        for (int j = i + 1; j < N; j++) {
            double c = A[j][i] / A[i][i];
            for (int k = i + 1; k < N; k++) A[j][k] -= c * A[i][k];
            b[j] -= c * b[i];
            A[j][i] = 0;
        }
    }
    for (int i = N - 1; i > -1; i--) {
        x[i] = b[i] / A[i][i];
        for (int j = i + 1; j < N; j++) x[i] -= A[i][j] / A[i][i] * x[j];
    }
}



/**
 * Метод прогонки
 */
void denisovrv::lab3()
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
    delete[] alpha;
    delete[] beta;

}



/**
 * Метод Холецкого
 */
void denisovrv::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void denisovrv::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void denisovrv::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void denisovrv::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void denisovrv::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void denisovrv::lab9()
{

}


std::string denisovrv::get_name()
{
  return "R.V. Denisov";
}

#include "akaykinsv.h"

/**
 * Введение в дисциплину
 */
void akaykinsv::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void akaykinsv::lab2()
{
    for (int column = 0, row = 0; column < N; column++, row++) {
        int max = row;
        for (int i = row + 1; i < N; i++) {
            if (abs(A[i][column]) > abs(A[max][column])) {
                max = i;
            }
        }

        if (max != row) {
            double swap[N];
            memcpy(swap, A[row], N);
            memcpy(A[row], A[max], N);
            memcpy(A[max], swap, N);
            double bSwap = b[max];
            b[max] = b [row];
            b[row] = bSwap;
        }
        double tmp = A[row][column];
        for (int i = column; i < N; i++) {
            A[row][i] /= tmp;
        }
        b[row] /= tmp;

        for (int i = row + 1; i < N; i++) {
            tmp = A[i][column];
            for (int j = column+1; j < N; j++) {
                A[i][j] -= A[row][j] * tmp;
            }
            b[i] -= b[row] * tmp;
        }
    } // получили треугольную матрицу

    for (int i = N-1, j = N - 1; i >= 0; i--, j--) {

        b[i] /= A[i][j];

        for (int k = i - 1; k >= 0; k--) {
            b[k] -= b[i] * A[k][j];
        }
        x[i] = b[i];
    }


}



/**
 * Метод прогонки
 */
void akaykinsv::lab3()
{

}



/**
 * Метод Холецкого
 */
void akaykinsv::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void akaykinsv::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void akaykinsv::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void akaykinsv::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void akaykinsv::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void akaykinsv::lab9()
{

}


std::string akaykinsv::get_name()
{
  return "S. V. Akaikin";
}

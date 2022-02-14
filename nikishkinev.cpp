#include "nikishkinev.h"

/**
 * Введение в дисциплину
 */
void nikishkinev::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void nikishkinev::lab2()
{
    int swaps[2 * N];
    for (int i = 0; i < N; i++) {
        int m = i;
        for (int j = i + 1; j < N; j++) {
            if (abs(A[m][i]) < abs(A[j][i])) {
                m = j;
            }
        }
        swap(A[i], A[m]);
        swap(b[i], b[m]);
        swaps[2 * i] = i;
        swaps[2 * i + 1] = m;
        for (int j = i + 1; j < N; j++) {
            A[i][j] /= A[i][i];
        }
        b[i] /= A[i][i];
        A[i][i] = 1;
        for (int j = i + 1; j < N; j++) {
            for (int k = i + 1; k < N; k++) {
                A[j][k] -= A[i][k] * A[j][i];
            }
            b[j] -= b[i] * A[j][i];
            A[j][i] = 0;
        }
    }

    for (int i = N - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            b[j] -= b[i] * A[j][i];
            A[j][i] = 0;
        }
    }

    for (int i = 2 * N - 1; i >= 0; i -= 2) {
        swap(A[swaps[i]], A[swaps[i - 1]]);
        swap(b[swaps[i]], b[swaps[i - 1]]);
    }

    for (int i = 0; i < N; i++) {
        x[i] = b[i];
    }
}



/**
 * Метод прогонки
 */
void nikishkinev::lab3()
{

}



/**
 * Метод Холецкого
 */
void nikishkinev::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void nikishkinev::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void nikishkinev::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void nikishkinev::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void nikishkinev::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void nikishkinev::lab9()
{

}


std::string nikishkinev::get_name()
{
  return "E.V. Nikishkin";
}

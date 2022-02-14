#include "ryabikinks.h"

/**
 * Введение в дисциплину
 */
void ryabikinks::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void ryabikinks::lab2()
{
    int i = 0;
    while (i < N) {
        int maxColumn = i;
        for (int j = i + 1; j < N; j++) {
            if (abs(A[j][i]) > abs(A[maxColumn][i])) {
                maxColumn = j;
            }
        }
        if (maxColumn != i) {
            swap(A[i], A[maxColumn]);
            swap(b[i], b[maxColumn]);
        }
        for (int j = i + 1; j < N; j++) {
            double ratio = A[j][i] / A[i][i];
            for (int k = i; k < N; k++) {
                A[j][k] -= (ratio * A[i][k]);
            }
            b[j] -= (ratio * b[i]);
        }
        i++;
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
void ryabikinks::lab3()
{

}



/**
 * Метод Холецкого
 */
void ryabikinks::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void ryabikinks::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void ryabikinks::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void ryabikinks::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void ryabikinks::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void ryabikinks::lab9()
{

}


std::string ryabikinks::get_name()
{
  return "K.S.Ryabikin";
}

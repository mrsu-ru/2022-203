#include "rodkinav.h"

/**
 * Введение в дисциплину
 */
void rodkinav::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void rodkinav::lab2()
{
    double Q = 0;

    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            Q = A[j][i] / A[i][i];

            for (int k = i; k < N; k++)
            {
                A[j][k] -= Q * A[i][k];
            }
            b[j] -= Q * b[i];
        }
    }

    for (int i = 0; i < N; i++)
    {
        x[i] = b[i];
    }

    for (int i = N - 1; i >= 0; i--)
    {
        for (int j = i + 1; j < N; j++)
        {
            x[i] -= A[i][j] * x[j];
        }

        x[i] /= A[i][i];
    }
}



/**
 * Метод прогонки
 */
void rodkinav::lab3()
{

}



/**
 * Метод Холецкого
 */
void rodkinav::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void rodkinav::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void rodkinav::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void rodkinav::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void rodkinav::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void rodkinav::lab9()
{

}


std::string rodkinav::get_name()
{
  return "A.V. Rodkin";
}

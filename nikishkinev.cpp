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
    double alpha[N];
    double beta[N];
    alpha[0] = A[0][1] / A[0][0];
    beta[0] = b[0] / A[0][0];

    for (int i = 1; i < N; i++) {
        alpha[i] = A[i][i + 1] / (A[i][i] - A[i][i - 1] * alpha[i - 1]);
        beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / (A[i][i] - A[i][i - 1] * alpha[i - 1]);
    }

    x[N - 1] = beta[N - 1];

    for (int i = N - 2; i >= 0; i--) {
        x[i] = beta[i] - alpha[i] * x[i + 1];
    }

}



/**
 * Метод Холецкого
 */
void nikishkinev::lab4()
{
    double S[N][N];
    double D[N];

    for (int i = 0; i < N; i++) {
        for (int l = 0; l < i; l++) {
            A[i][i] -= S[l][i] * S[l][i] * D[l];
        }
        D[i] = A[i][i] >= 0 ? 1 : -1;
        S[i][i] = sqrt(D[i] * A[i][i]);

        for (int j = i + 1; j < N; j++) {
            for (int l = 0; l < j; l++) {
                A[i][j] -= S[l][i] * D[l] * S[l][j];
            }
            S[i][j] = A[i][j] / (S[i][i] * D[i]);
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
            b[i] -= S[j][i] * b[j] * D[j];
        }
        b[i] /= S[i][i] * D[i];
    }

    for (int i = N - 1; i >= 0; i--) {
        for (int j = i + 1; j < N; j++) {
            b[i] -= S[i][j] * x[j];
        }
        x[i] = b[i] / S[i][i];
    }


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

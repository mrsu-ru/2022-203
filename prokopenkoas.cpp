#include "prokopenkoas.h"

/**
 * Введение в дисциплину
 */
void prokopenkoas::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void prokopenkoas::lab2()
{
    int max;
    double sum = 0;
    //прямой ход
    for (int i = 0; i < N; i++)
    {
        max = i;
        double* vspom; double vspom_1;

        for (int j = i + 1; j < N; j++)
        {
            if (fabs(A[j][i]) > fabs(A[max][i]))
            {
                max = j;
            }
        }
        if (max != i)
        {
            vspom = A[i];
            A[i] = A[max];
            A[max] = vspom;

            vspom_1 = b[i];
            b[i] = b[max];
            b[max] = vspom_1;
        }

        for (int j = N - 1; j > i; j--)
        {
            A[i][j] /= A[i][i];
        }
        b[i] /= A[i][i];
        A[i][i] = 1;

        for (int j = i + 1; j < N; j++)
        {
            for (int k = N - 1; k > i; k--)
            {
                A[j][k] = A[j][k] - A[j][i] * A[i][k];
            }
            b[j] = b[j] - A[j][i] * b[i];
            A[j][i] = 0;
        }
    }

    //обратный ход
    x[N - 1] = b[N - 1];
    for (int i = N - 2; i > -1; i--)
    {
        for (int j = i + 1; j < N; j++)
        {
            sum += x[j] * A[i][j];
        }
        x[i] = b[i] - sum;
        sum = 0;
    }
}



/**
 * Метод прогонки
 */
void prokopenkoas::lab3()
{
    double* alfa = new double[N];
    double* beta = new double[N];
    alfa[0] = A[0][1] / -A[0][0];
    beta[0] = b[0] / A[0][0];
    for (int i = 1; i < N; i++) {
        alfa[i] = A[i][i + 1] / (-A[i][i] - A[i][i - 1] * alfa[i - 1]);
        beta[i] = (A[i][i - 1] * beta[i - 1] - b[i]) / (-A[i][i] - A[i][i - 1] * alfa[i - 1]);
    }
    x[N - 1] = beta[N - 1];
    for (int i = N - 2; i >= 0; i--) {
        x[i] = alfa[i] * x[i + 1] + beta[i];
    }
}



/**
 * Метод Холецкого
 */
void prokopenkoas::lab4()
{
    double** S = new double* [N];
    for (int i = 0; i < N; i++) {
        S[i] = new double[N];
        for (int j = 0; j < N; j++) S[i][j] = 0;
    }

    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if (j == i) {
                if (i == 0) {
                    S[i][i] = sqrt(A[0][0]);
                }
                else {
                    for (int k = 0; k < i; k++) A[i][i] -= S[k][i] * S[k][i];
                    S[i][i] = sqrt(A[i][i]);
                }
            }
            else {
                for (int k = 0; k < j; k++) A[i][j] -= S[k][i] * S[k][j];
                S[i][j] = A[i][j] / S[i][i];
            }
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) b[i] -= S[j][i] * b[j];
        b[i] /= S[i][i];
    }

    for (int i = N - 1; i >= 0; i--) {
        for (int k = i + 1; k < N; k++) b[i] -= S[i][k] * x[k];
        x[i] = b[i] / S[i][i];
    }

    for (int i = 0; i < N; i++) delete[] S[i];
    delete[] S;


}



/**
 * Метод Якоби или Зейделя
 */
void prokopenkoas::lab5()
{
    double norm, sum;
    double eps = 1.e-15;
    for (int i = 0; i < N; ++i) x[i] = b[i];

    do {
        norm = 0;
        for (int i = 0; i < N; i++) {
            sum = 0;
            for (int j = 0; j < N; j++) {
                if (i != j) sum += A[i][j] * x[j];
            }
            sum = (b[i] - sum) / A[i][i];
            if (norm < (fabs(sum - x[i]))) {
                norm = (fabs(sum - x[i]));
            }
            x[i] = sum;
        }
    } while (norm >= eps);
}


double* MulVecToMatrix(int N, double* A[], double b[]) {
    double* temp = new double[N];
    for (int i = 0; i < N; i++) {
        temp[i] = 0;
        for (int j = 0; j < N; j++) {
            temp[i] += A[i][j] * b[i];
        }
    }
    return temp;
}

double ScalarMul(int N, double temp[], double r[]) {
    double k = 0;
    for (int i = 0; i < N; i++) {
        k += temp[i] * r[i];
    }
    return k;

}



/**
 * Метод минимальных невязок
 */
void prokopenkoas::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void prokopenkoas::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void prokopenkoas::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void prokopenkoas::lab9()
{

}


std::string prokopenkoas::get_name()
{
  return "A.S. Prokopenko";
}

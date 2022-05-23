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
    double S[N][N];
    int D[N][N];
    for (int i = 0; i < N; i++) {

        double Sum = 0;
        for (int l = 0; l <= i - 1; l++) {
            Sum += S[l][i] * S[l][i] * D[l][l];
        }

        int sign = (A[i][i] - Sum) < 0;
        D[i][i] = (int)pow(-1, sign);

        S[i][i] = sqrt(abs(A[i][i] - Sum));
        for (int j = 0; j < i; j++) S[i][j] = 0.0;
        for (int j = i + 1; j < N; j++) {
            double Sum = 0;
            for (int l = 0; l <= i - 1; l++) Sum += S[l][j] * S[l][i] * D[l][l];
            S[i][j] = (A[i][j] - Sum) / S[i][i] * D[i][i];
        }

        // обратно

        double y[N];
        y[0] = b[0] / S[0][0];
        for (int i = 1; i < N; i++) {
            double Sum = 0;
            for (int j = 0; j <= i - 1; j++) Sum += S[j][i] * y[j];
            y[i] = b[i] - Sum;
            y[i] /= S[i][i];
        }
        x[N - 1] = y[N - 1] / S[N - 1][N - 1];
        for (int i = N - 2; i >= 0; i--) {
            double Sum = 0;
            for (int k = i + 1; k <= N - 1; k++) Sum += S[i][k] * x[k];
            x[i] = y[i] - Sum;
            x[i] /= S[i][i];
        }

    }

}



/**
 * Метод Якоби или Зейделя
 */
void prokopenkoas::lab5()
{
    double* toka = new double[N];
    double eps = 1e-20;
    double norm;
    for (int i = 0; i < N; i++) x[i] = b[i] / A[i][i];

    do {
        for (int i = 0; i < N; i++) {

            for (int j = 0; j < N; j++) toka[j] = x[j];

            double lowSum = 0, uppSum = 0;

            for (int j = 0; j < i; j++) lowSum += A[i][j] * toka[j];
            for (int j = i + 1; j < N; j++) uppSum += A[i][j] * toka[j];

            x[i] = 1 / A[i][i] * (b[i] - lowSum - uppSum);

            if (i == 0) norm = abs(x[i] - toka[i]);
            if (abs(x[i] - toka[i]) > norm) norm = abs(x[i] - toka[i]);

        }

    } while (norm >= eps);
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

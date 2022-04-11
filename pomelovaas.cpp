#include "pomelovaas.h"

/**
 * Введение в дисциплину
 */
void pomelovaas::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void pomelovaas::lab2()
{
    double* temp = b;
    for (int i = 0, j = 0; (i < N - 1) && (j < N - 1); i++, j++)
    {
        int cur = i;
        double max = A[j][j];
        for (int num = j; num < N; num++)
            if (abs(max) < abs(A[num][j])) max = A[num][j];

        while (A[i][j] != max)
            i++;

        if (cur != i)
        {
            swap(temp[cur], temp[i]);
            swap(A[cur], A[i]);
            i = cur;
        }

        double t = 1/A[i][j];
        temp[i] *= t;
        for (int cj = j; cj < N; cj++)
            A[i][cj] *= t;

        for (int ci = i + 1; ci < N; ci++)
        {
            temp[ci] += temp[i] * (-A[ci][j]);
            t = -A[ci][j];
            for (int cj = j; cj < N; cj++)
            {
                A[ci][cj] += A[i][cj] * t;
            }
        }
    }

    for (int i = N - 1, j = N - 1; (i > 0) && (j > 0); i--, j--)
    {
        ///int cur = i;
        temp[i] *= 1 / A[i][j];
        A[i][j] *= 1 / A[i][j];
        for (int ci = i - 1; ci >= 0; ci--)
        {
            temp[ci] += temp[i] * (-A[ci][j]);
            A[ci][j] = 0;
        }
    }

    x = temp;
}

/**
 * Метод прогонки
 */
void pomelovaas::lab3()
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
void pomelovaas::lab4()
{
    double **S = new double* [N];
    for (int i = 0; i < N; ++i) S[i] = new double[N];
    S[0][0] = sqrt( fabs(A[0][0]) );
    double n; //
    double *y = new double[N];
    double *D = new double[N];
    A[0][0] > 0 ?
            D[0] = 1
                :
            D[0] = -1;


    for (int i = 1; i < N; ++i) S[0][i] = A[0][i] / ( D[0]*S[0][0] );

    for (int i = 1; i < N; ++i){
        n = 0;
        for (int j=0; j<i; j++) n += D[j] * S[j][i] * S[j][i];

        A[i][i] - n >= 0 ?
                D[i] = 1
                         :
                D[i] = -1;

        S[i][i] = sqrt( D[i]*(A[i][i] - n) );

        for (int j= i+1; j < N; ++j) {
            double l = 0;
            for (int k = 0; k < j; ++k) l += D[k] * S[k][i] * S[k][j];

            S[i][j] =
                    ( A[i][j]-l )
                    /
                    ( D[i]*S[i][i] );
        }

    }


    for (int i = 0; i < N; ++i){
        n=0;
        if(i == 0){
            y[0] = b[0] / S[0][0];
            ++i;
        }

        for (int j=0; j < i; ++j){
            n += y[j]*S[j][i];
        }
        y[i] = (b[i] - n) / S[i][i];
    }

    x[N-1] = y[N-1] / (D[N-1] * S[N-1][N-1]);

    for (int i = N-2; i >= 0; --i){
        n = 0;
        for (int j = i+1; j < N; ++j) n += x[j] * D[j] * S[i][j];

        x[i]=
                (y[i] - n)
                /
                ( D[i]*S[i][i] );
    }
}




/**
 * Метод Якоби или Зейделя
 */
void pomelovaas::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void pomelovaas::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void pomelovaas::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void pomelovaas::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void pomelovaas::lab9()
{

}


std::string pomelovaas::get_name()
{
  return "A.S.Pomelova";
}

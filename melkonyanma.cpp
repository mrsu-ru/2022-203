#include "melkonyanma.h"
#include <iostream>
#include <math.h>

#define epsilon 1e-20

/**
 * Введение в дисциплину
 */
void melkonyanma::lab1()
{
    cout << "Hello World!" << endl;
}

/**
 * Метод Гаусса с выбором главного элемента
 */
void melkonyanma::lab2()
{
    method_Gauss(A, b, N);

    for (int i = 0; i < N; ++i)
    {
        x[i] = b[i];
    }
}

/**
 * Метод прогонки
 */
void melkonyanma::lab3()
{
    double *alpha = new double[N];
    double *beta = new double[N];
    double kappa;

    init_matrix(alpha, N);
    init_matrix(beta, N);

    kappa = A[0][0];
    alpha[0] = -A[0][1] / kappa;
    beta[0] = b[0] / kappa;

    for (int i = 1; i < N; ++i)
    {
        kappa = A[i][i] + A[i][i - 1] * alpha[i - 1];

        if (i != N - 1)
        {
            alpha[i] = -A[i][i + 1] / kappa;
        }

        beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / kappa;
    }

    x[N - 1] = beta[N - 1];

    for (int i = N - 2; i >= 0; --i)
    {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }

    delete[] alpha;
    delete[] beta;
}

/**
 * Метод Холецкого
 */
void melkonyanma::lab4()
{
    double *y = new double[N];
    double **L = new double *[N];
    double **Lt = new double *[N];
    double summ;

    create_dynamic_array(L, N);
    create_dynamic_array(Lt, N);

    init_matrix(L, N);
    init_matrix(Lt, N);
    init_matrix(y, N);

    L[0][0] = sqrt(A[0][0]);
    for (int i = 1; i < N; ++i)
    {
        L[i][0] = A[i][0] / L[0][0];
    }

    for (int i = 1; i < N; ++i)
    {
        summ = 0.;

        for (int p = 0; p < i; ++p)
        {
            summ += L[i][p] * L[i][p];
        }
        L[i][i] = sqrt(A[i][i] - summ);

        summ = 0.;

        for (int j = i + 1; j < N; ++j)
        {
            for (int p = 0; p < i; ++p)
            {
                summ += L[i][p] * L[j][p];
            }

            L[j][i] = (A[j][i] - summ) / L[i][i];
        }
    }

    MatrixTrans(L, Lt, N);
    method_Gauss(L, b, N);
    for (int i = 0; i < N; ++i)
    {
        y[i] = b[i];
    }
    method_Gauss(Lt, y, N);
    for (int i = 0; i < N; ++i)
    {
        x[i] = y[i];
    }

    delete_dynamic_array(L, N);
    delete_dynamic_array(Lt, N);

    delete[] y;
}

/**
 * Метод Якоби или Зейделя
 */
void melkonyanma::lab5()
{
    double xk[N];
    double norm;
    double LU[N][N];

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (i != j)
            {
                LU[i][j] = A[i][j];
            }

            if (i == j)
            {
                LU[i][j] = 0.;
            }
        }
    }

    for (int i = 0; i < N; ++i)
        x[i] = b[i] / A[i][i];

    do
    {
        for (int i = 0; i < N; ++i)
        {
            for (int i = 0; i < N; ++i)
                xk[i] = x[i];

            double summ = 0.;

            for (int j = 0; j < N; ++j)
            {
                summ += LU[i][j] * xk[j];
            }

            x[i] = (1 / A[i][i]) * (b[i] - summ);

            if (i == 0)
            {
                norm = fabs(x[i] - xk[i]);
            }

            if (fabs(x[i] - xk[i]) > norm)
            {
                norm = fabs(x[i] - xk[i]);
            }
        }
    } while (sqrt(norm) >= epsilon);
}

/**
 * Метод минимальных невязок
 */
void melkonyanma::lab6()
{
    double norm = 1e-1;
    double rk[N];
    init_matrix(rk, N);

    for (int i = 0; i < N; ++i)
    {
        x[i] = 1e-2;
    }

    do
    {
        double *Ax = MatrixMultOnVect(A, x, N);

        for (int i = 0; i < N; ++i)
        {
            rk[i] = Ax[i] - b[i];
        }

        double *Ar = MatrixMultOnVect(A, rk, N);
        double tau = 0;

        tau = (scalar(Ar, rk, N) / scalar(Ar, Ar, N));

        for (int i = 0; i < N; i++)
        {
            double prevX = x[i];
            x[i] = prevX - tau * rk[i];

            if (fabs(x[i] - prevX) < norm)
            {
                norm = fabs(x[i] - prevX);
            }
        }

        delete[] Ar;
        delete[] Ax;

    } while (sqrt(norm) >= epsilon);
}

/**
 * Метод сопряженных градиентов
 */
void melkonyanma::lab7() {}

/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void melkonyanma::lab8()
{
    double eps = 1e-1;
    double max_element;
    double c, s;
    int posMaxElementI, posMaxElementJ;

    do
    {
        max_element = 0.;

        for (int i = 0; i < N - 1; ++i)
        {
            for (int j = i + 1; j < N; ++j)
            {
                if (abs(A[i][j]) > max_element)
                {
                    max_element = abs(A[i][j]);
                    posMaxElementI = i;
                    posMaxElementJ = j;
                }
            }
        }

        double aii, ajj;
        double teta = 0.;

        aii = A[posMaxElementI][posMaxElementI];
        ajj = A[posMaxElementJ][posMaxElementJ];

        if ((aii - ajj) < eps)
        {
            teta = atan(1);
        }
        else
        {
            teta = atan((2 * max_element) / (aii - ajj)) / 2;
        }

        c = cos(teta);
        s = sin(teta);

        double **H = new double *[N];
        create_dynamic_array(H, N);
        double **transH = new double *[N];
        create_dynamic_array(transH, N);

        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                if (i == j)
                    H[i][j] = 1;
                else
                    H[i][j] = 0;
            }
        }

        H[posMaxElementI][posMaxElementI] = H[posMaxElementJ][posMaxElementJ] = c;
        H[posMaxElementI][posMaxElementJ] = -s;
        H[posMaxElementJ][posMaxElementI] = s;

        MatrixTrans(H, transH, N);

        A = MatrixMultOnMatrix(MatrixMultOnMatrix(transH, A, N), H, N);

        delete_dynamic_array(H, N);
        delete_dynamic_array(transH, N);
        delete[] H;
        delete[] transH;

    } while (max_element >= eps);

    for (int i = 0; i < N; ++i)
    {
        printf("Lambda(%d) = %.4f\n", i + 1, A[i][i]);
        x[i] = A[i][i];
    }
}

/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void melkonyanma::lab9()
{
    double AbsMaxEigenvalue;
    double Eigenvalue = 0.;
    double eps = 1e-8;
    double *xPrev = new double[N];
    double *xNew = new double[N];

    for (int i = 0; i < N; ++i)
    {
        xPrev[i] = b[i];
    }

    do
    {
        AbsMaxEigenvalue = Eigenvalue;

        xNew = MatrixMultOnVect(A, xPrev, N);

        double s1 = 0, s2 = 0;
        for (int i = 0; i < N; ++i)
        {
            s1 += xNew[i];
            s2 += xPrev[i];
        }

        Eigenvalue = s1 / s2;

        for (int i = 0; i < N; ++i)
        {
            xPrev[i] = xNew[i] / scalar(xNew, xNew, N);
        }

    } while (fabs(AbsMaxEigenvalue - Eigenvalue) >= eps);

    delete[] xNew;
    delete[] xPrev;

    cout << "Absolute max eigenvalue = " << AbsMaxEigenvalue << endl;
}

std::string melkonyanma::get_name()
{
    return "M.A. Melkonyan";
}
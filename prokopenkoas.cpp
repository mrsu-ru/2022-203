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

        double subSum = 0;
        for (int l = 0; l <= i - 1; l++) {
            subSum += S[l][i] * S[l][i] * D[l][l];
        }

        int sign = (A[i][i] - subSum) < 0;
        D[i][i] = (int)pow(-1, sign);

        S[i][i] = sqrt(abs(A[i][i] - subSum));
        for (int j = 0; j < i; j++) S[i][j] = 0.0;
        for (int j = i + 1; j < N; j++) {
            double subSum = 0;
            for (int l = 0; l <= i - 1; l++) subSum += S[l][j] * S[l][i] * D[l][l];
            S[i][j] = (A[i][j] - subSum) / S[i][i] * D[i][i];
        }

        //Обратный ход

        double y[N];
        y[0] = b[0] / S[0][0];
        for (int i = 1; i < N; i++) {
            double subSum = 0;
            for (int j = 0; j <= i - 1; j++) subSum += S[j][i] * y[j];
            y[i] = b[i] - subSum;
            y[i] /= S[i][i];
        }
        x[N - 1] = y[N - 1] / S[N - 1][N - 1];
        for (int i = N - 2; i >= 0; i--) {
            double subSum = 0;
            for (int k = i + 1; k <= N - 1; k++) subSum += S[i][k] * x[k];
            x[i] = y[i] - subSum;
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


double* MulVecToMatrica(int N, double* A[], double b[]) {
    double* temp = new double[N];
    for (int i = 0; i < N; i++) {
        temp[i] = 0;
        for (int j = 0; j < N; j++) {
            temp[i] += A[i][j] * b[i];
        }
    }
    return temp;
}

double ScalarMullet(int N, double temp[], double r[]) {
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
    double* r = new double[N];
    double* xk = new double[N];
    double eps = 1.e-17;
    double* x1 = b;
    double t;
    double maxDelta;
    do {
        double* temp = MulVecToMatrica(N, A, x1);

        for (int i = 0; i < N; i++) {
            r[i] = temp[i] - b[i];
        }
        double* Ar = MulVecToMatrica(N, A, r);

        double Scalar1, Scalar2;
        Scalar1 = ScalarMullet(N, Ar, r);
        Scalar2 = ScalarMullet(N, Ar, Ar);
        t = Scalar1 / Scalar2;

        for (int i = 0; i < N; i++) {
            xk[i] = x1[i] - t * r[i];
        }
        maxDelta = abs(xk[0] - x1[0]);

        for (int i = 1; i < N; i++) {
            double delta = abs(xk[i] - x1[i]);
            if (delta > maxDelta) {
                maxDelta = delta;
            }
        }
        x1 = xk;

    } while (maxDelta > eps);

    for (int i = 0; i < N; i++) {
        x[i] = x1[i];
    }
    delete[] r;
    delete[] xk;
}



/**
 * Метод сопряженных градиентов
 */
void prokopenkoas::lab7()
{

}


double** TransposeMatrica(double**& m, int n) {
    double** temp = new double* [n];
    for (int i = 0; i < n; ++i) {
        temp[i] = new double[n];
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            temp[i][j] = m[j][i];
        }
    }

    return temp;
}


double** MulMatricaa(double**& m1, double**& m2, int n) {
    double** temp = new double* [n];
    for (int i = 0; i < n; ++i) {
        temp[i] = new double[n];
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            temp[i][j] = 0;
            for (int k = 0; k < n; ++k) {
                temp[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }

    return temp;
}

/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void prokopenkoas::lab8()
{
    double eps = 1.e-1;
    double MaxEl;

    do {
        MaxEl = 0;
        int maxi = 0, maxj = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i < j) {
                    if (MaxEl < fabs(A[i][j])) {
                        MaxEl = fabs(A[i][j]);
                        maxi = i;
                        maxj = j;
                    }
                }
            }
        }

        double Phi = 0.5 * atan(2 * A[maxi][maxj] / (A[maxi][maxi] - A[maxj][maxj]));

        double** m_H = new double* [N];
        for (int i = 0; i < N; ++i) {
            m_H[i] = new double[N];
        }

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i == j) {
                    m_H[i][j] = 1;
                }
                else {
                    m_H[i][j] = 0;
                }
            }
        }

        m_H[maxi][maxj] = -sin(Phi);
        m_H[maxj][maxi] = sin(Phi);

        m_H[maxi][maxi] = cos(Phi);
        m_H[maxj][maxj] = cos(Phi);

        double** t_H = TransposeMatrica(m_H, N);

        double** t_A = MulMatricaa(t_H, A, N);

        A = MulMatricaa(t_A, m_H, N);
    } while (MaxEl >= eps);

    for (int i = 0; i < N; ++i) {
        cout << "lambda  = " << A[i][i] << endl;
    }
}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void prokopenkoas::lab9()
{
  
        const double eps = 1e-3;

        double* yPrev = new double[N];
        for (int i = 0; i < N; i++)
            yPrev[i] = 1;

        double* y = new double[N];

        int iter;
        double delta = 1e9;
        double lambdaMax = 0;

        for (iter = 0; delta > eps; iter++) {
            for (int i = 0; i < N; i++) {
                y[i] = 0;
                for (int j = 0; j < N; j++) {
                    y[i] += A[i][j] * yPrev[j];
                }
            }

            double sum1 = 0, sum2 = 0;
            for (int i = 0; i < N; ++i) {
                sum1 += y[i];
                sum2 += yPrev[i];
            }
            double lambda = sum1 / sum2;
            delta = fabs(lambda - lambdaMax);
            lambdaMax = lambda;

            for (int i = 0; i < N; i++) {
                yPrev[i] = y[i];
            }
            double r = 0;
            for (int i = 0; i < N; i++) {
                r += y[i] * y[i];
            }
            r = sqrt(r);
            for (int i = 0; i < N; ++i) {
                y[i] /= r;
            }
        }
        printf("Max Eigen value -- %.4f\n", lambdaMax);

        delete[]yPrev;
        delete[]y;
    
}


std::string prokopenkoas::get_name()
{
  return "A.S. Prokopenko";
}

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
    double alpha[N - 1], beta[N], y[N];
    //Прямая прогонка
    y[0] = A[0][0];
    beta[0] = b[0] / y[0];
    alpha[0] = -A[0][1] / y[0];
    for (int i = 1; i < N; i++) {
        y[i] = A[i][i] + A[i][i-1]*alpha[i-1];
        beta[i] = (b[i] - A[i][i-1]*beta[i-1])/y[i];
        if (i != N - 1) {
            alpha[i] = -A[i][i+1]/y[i];
        }
    }

    //Обратная прогонка
    x[N-1] = beta[N-1];
    for (int i = N-2; i >= 0; i--) {
        x[i] = alpha[i] * x[i+1] + beta[i];
    }
}


/**
 * Метод Холецкого
 */
void akaykinsv::lab4()
{
    double S[N][N];
    int D[N][N];
    for (int i = 0; i < N; i++) {

        double subSum = 0;
        for (int l = 0; l <= i - 1; l++) {
            subSum += S[l][i] * S[l][i] * D[l][l];
        }

        int sign = (A[i][i] - subSum) < 0;
        D[i][i] =  (int)pow(-1, sign);

        S[i][i] = sqrt(abs(A[i][i] - subSum));
        for (int j = 0; j < i; j++) S[i][j] = 0.0;
        for (int j = i+1; j < N; j++) {
            double subSum = 0;
            for (int l = 0; l <= i - 1; l++) subSum += S[l][j] * S[l][i] * D[l][l];
            S[i][j] = (A[i][j] - subSum)/S[i][i]*D[i][i];
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
        for (int i = N - 2; i >=0; i--) {
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
void akaykinsv::lab5()
{
    double *xk = new double[N];
    double eps = 1e-20;
    double norm;
    //Начальное приближение
    for (int i = 0; i < N; i++) x[i] = b[i]/A[i][i];

    do {
        for (int i = 0; i < N; i++) {

            for (int j = 0; j < N; j++) xk[j] = x[j];

            double lowerSum = 0, upperSum = 0;

            for (int j = 0; j < i; j++) lowerSum += A[i][j] * xk[j];
            for (int j = i + 1; j < N; j++) upperSum += A[i][j] * xk[j];

            x[i] = 1/A[i][i] * (b[i] - lowerSum - upperSum);

            if ( i == 0 ) norm = abs(x[i] - xk[i]);
            if (abs(x[i] - xk[i]) > norm) norm = abs(x[i] - xk[i]);

        }

    } while (norm >= eps);
}



/**
 * Метод минимальных невязок
 */
 // Вспомогательный метод произведения матриц на вектор
double *MVcomposition (double **A, double *x, int N) {
    double *result = new double[N];
    for (int i=0; i < N; i++) {
        result[i] = 0;
        for (int j = 0; j < N; j++) {
            result[i] += A[i][j] * x[j];
        }
    }
    return result;
}
//Вспомогательный метод вычисления скалярного произведения двух векторов
double scalar(double *x, double *y, int N) {
    double scal = 0;
    for (int i = 0; i < N; i++) scal+= x[i] * y[i];
    return scal;
}
void akaykinsv::lab6()
{
    double eps = 1e-19;
    double *rk = new double[N];
    double T = 0;
    for (int i = 0; i < N; i++) x[i] = 0.1;

    double norm;
    do {
        double *Axk = MVcomposition(A, x, N);

        for (int i = 0; i < N; i++) rk[i] = Axk[i] - b[i];

        double *Ark = MVcomposition(A, rk, N);

        T = scalar(Ark, rk, N) / scalar(Ark, Ark, N);
        norm = 0;
        for (int i = 0; i < N; i++) {
            double check = x[i];
            x[i] = x[i] - T * rk[i];
            norm += (x[i] - check) * (x[i] - check);
        }

    } while (sqrt(norm) > eps);


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
    double eps = 1e-20;
    double** B = new double* [N];
    for (int i = 0; i < N; i++) {
        B[i] = new double[N];
    }

    while (true) {
        double norm = 0;
        int imax = 0;
        int jmax = 1;
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                if (abs(A[i][j]) > abs(A[imax][jmax])) {
                    imax = i;
                    jmax = j;
                }
                norm += A[i][j] * A[i][j];
            }
        }

        if (sqrt(norm) < eps) {
            break;
        }

        double fi = 0.5 * atan(2 * A[imax][jmax] / (A[imax][imax] - A[jmax][jmax]));

        for (int i = 0; i < N; i++) {
            B[i][imax] = A[i][imax] * cos(fi) + A[i][jmax] * sin(fi);
            B[i][jmax] = A[i][jmax] * cos(fi) - A[i][imax] * sin(fi);
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (j != imax && j != jmax) {
                    B[i][j] = A[i][j];
                }
            }
        }

        for (int j = 0; j < N; j++) {
            A[imax][j] = B[imax][j] * cos(fi) + B[jmax][j] * sin(fi);
            A[jmax][j] = B[jmax][j] * cos(fi) - B[imax][j] * sin(fi);
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i != imax && i != jmax) {
                    A[i][j] = B[i][j];
                }
            }
        }
    }

    for (int i = 0; i < N; i++) {
        x[i] = A[i][i];
    }
}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void akaykinsv::lab9()
{
    double eps = 1e-3;
    double *yk = new double[N];
    double *yk1 = new double[N];
    int i,j;
    for (int i = 0; i < N; i++) yk[i] = 1;
    double err = 0;
    double lambda1 = 0;
    double lambda2 = 0;
    do{
        for (i = 0; i < N; i++)
            for (j = 0; j < N; j++)
                yk1[i] += A[i][j] * yk[j];
        for (i = 0; i < N; i++)
            if(fabs(yk[i]) > eps && fabs(yk1[i]) > eps){
                lambda1 = yk1[i]/yk[i];
                break;
            }
        err = fabs(lambda1 - lambda2);
        lambda2 = lambda1;
        for (i = 0; i < N; i++) {
            yk[i] = yk1[i];
            yk1[i] = 0;
        }
    }
    while(err > eps);
    cout<<"max lambda = "<<lambda2<<endl;
}



std::string akaykinsv::get_name()
{
  return "S. V. Akaikin";
}

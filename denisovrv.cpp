#include "denisovrv.h"

/**
 * Введение в дисциплину
 */
void denisovrv::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void denisovrv::lab2()
{
    for (int i = 0; i < N; i++) {
        int max = 0;
        for (int j = 0; j < N; j++) if (abs(A[j][i]) > abs(A[max][i])) max = j;
        swap(A[i], A[max]);
        swap(b[i], b[max]);
        for (int j = i + 1; j < N; j++) {
            double c = A[j][i] / A[i][i];
            for (int k = i + 1; k < N; k++) A[j][k] -= c * A[i][k];
            b[j] -= c * b[i];
            A[j][i] = 0;
        }
    }
    for (int i = N - 1; i > -1; i--) {
        x[i] = b[i] / A[i][i];
        for (int j = i + 1; j < N; j++) x[i] -= A[i][j] / A[i][i] * x[j];
    }
}



/**
 * Метод прогонки
 */
void denisovrv::lab3()
{
    double* alpha = new double[N];
    double* beta = new double[N];
    double y;

    y = A[0][0];
    alpha[0] = -A[0][1] / y;
    beta[0] = b[0] / y;

    for (int i = 1; i < N; i++) {
        y = A[i][i] + A[i][i - 1] * alpha[i - 1];
        alpha[i] = -A[i][i + 1] / y;
        beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / y;
    }
    x[N - 1] = beta[N - 1];

    for (int i = N - 2; i >= 0; i--) {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }
    delete[] alpha;
    delete[] beta;

}



/**
 * Метод Холецкого
 */
void denisovrv::lab4()
{
    double** L = new double* [N];
    for (int i = 0; i < N; i++) {
        L[i] = new double[N];
        for (int j = 0; j < N; j++) L[i][j] = 0;
    }
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if (j == i) {
                for (int l = 0; l < i; l++) A[i][i] -= L[l][i] * L[l][i];
                L[i][i] = sqrt(A[i][i]);
            }
            else {
                for (int l = 0; l < j; l++) A[i][j] -= L[l][i] * L[l][j];
                L[i][j] = A[i][j] / L[i][i];
            }
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) b[i] -= L[j][i] * b[j];
        b[i] /= L[i][i];
    }
    for (int i = N - 1; i >= 0; i--) {
        for (int j = i + 1; j < N; j++) b[i] -= L[i][j] * x[j];
        x[i] = b[i] / L[i][i];
    }
    for (int i = 0; i < N; i++) delete[] L[i];
    delete[] L;
}



/**
 * Метод Якоби или Зейделя
 */
void denisovrv::lab5()
{
    double norm, sum, appr;
    double eps = 1.e-15;
    for (int i = 0; i < N; ++i) x[i] = b[i];
    do {
        norm = 0;
        for (int i = 0; i < N; i++) {
            sum = 0;
            for (int j = 0; j < N; j++) if (i != j) sum += A[i][j] * x[j];
            appr = (b[i] - sum) / A[i][i];
            if (norm < (abs(appr - x[i]))) norm = (abs(appr - x[i]));
            x[i] = appr;
        }
    } while (norm >= eps);
}



/**
 * Метод минимальных невязок
 */
void denisovrv::lab6()
{
    double* x1 = new double[N];
    double* r = new double[N];
    double* Ar = new double[N];
    double eps = 1.e-15;
    double Ar2 = 0;
    double Arr = 0;
    double t = 0;
    double norma = 0;
    int i, j;

    for (i = 0; i < N; i++) {
        x[i] = b[i];
        x1[i] = 0;
        Ar[i] = 0;
    }

    for (i = 0; i < N; i++) {
        r[i] = b[i];
        for (j = 0; j < N; j++) r[i] -= A[i][j] * x[j];
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) Ar[i] += A[i][j] * r[j];
    }

    for (i = 0; i < N; i++) {
        Ar2 += Ar[i] * Ar[i];
        Arr += Ar[i] * r[i];
    }

    t = -Arr / Ar2;

    for (i = 0; i < N; i++) x1[i] = x[i] - t * r[i];

    norma = abs(x1[0] - x[0]);
    for (i = 1; i < N; i++) {
        if ((abs(x1[i] - x[i])) > norma) norma = abs(x1[i] - x[i]);
    }

    while (norma > eps) {
        for (i = 0; i < N; i++) {
            x[i] = x1[i];
            x1[i] = 0;
            Ar[i] = 0;
            r[i] = 0;
        }

        Ar2 = 0;
        Arr = 0;
        t = 0;
        norma = 0;

        for (i = 0; i < N; i++) {
            r[i] = b[i];
            for (j = 0; j < N; j++) r[i] -= A[i][j] * x[j];
        }

        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) Ar[i] += A[i][j] * r[j];
        }

        for (i = 0; i < N; i++) {
            Ar2 += Ar[i] * Ar[i];
            Arr += Ar[i] * r[i];
        }

        t = -Arr / Ar2;

        for (i = 0; i < N; i++) x1[i] = x[i] - t * r[i];

        norma = abs(x1[0] - x[0]);
        for (i = 1; i < N; i++) {
            if ((abs(x1[i] - x[i])) > norma) norma = abs(x1[i] - x[i]);
        }
    }

    for (i = 0; i < N; i++) x[i] = x1[i];

    delete[]x1;
    delete[]r;
    delete[]Ar;
}



/**
 * Метод сопряженных градиентов
 */
void denisovrv::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void denisovrv::lab8()
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
void denisovrv::lab9()
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            if (A[i][j] != A[j][i]) return;
    }
    double l, maxl = 0;
    bool flag = true;
    double* y = new double[N];
    double* y_prev = new double[N];
    for (int i = 0; i < N; i++) y_prev[i] = 1;
    while (flag) {
        flag = false;
        for (int i = 0; i < N; i++) {
            y[i] = 0;
            for (int j = 0; j < N; j++) {
                y[i] += A[i][j] * y_prev[j];
            }
        }
        l = 0;
        for (int i = 0; i < N; i++) {
            if (fabs(y[i]) > eps && fabs(y_prev[i]) > eps) {
                l = y[i] / y_prev[i];
                break;
            }
        }
        if (fabs(l - maxl) > eps) flag = true;
        maxl = l;
        for (int i = 0; i < N; i++) y_prev[i] = y[i];
    }
    cout << "The largest lambda: " << maxl << endl;
    delete[] y;
    delete[] y_prev;
}


std::string denisovrv::get_name()
{
  return "R.V. Denisov";
}

#include "venediktovayap.h"

/**
 * Введение в дисциплину
 */
void venediktovayap::lab1() {
    cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void venediktovayap::lab2() {

    for (int i = 0; i < N; i++) {
        int maxLine = i;
        for (int j = i + 1; j < N; j++) {
            if (abs(A[j][i]) > abs(A[maxLine][i])) {
                maxLine = j;
            }
        }
        if (i != maxLine) {
            swap(A[maxLine], A[i]);
            swap(b[maxLine], b[i]);
        }
        for (int j = i + 1; j < N; j++) {
            double c = A[j][i] / A[i][i];
            for (int k = i; k < N; k++) {
                A[j][k] -= (c * A[i][k]);
            }
            b[j] -= (c * b[i]);
        }
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
void venediktovayap::lab3() {
    double alpha[N];
    double beta[N];
    double c;

    c = A[0][0];
    alpha[0] = -A[0][1] / c;
    beta[0] = b[0] / c;
    for (int i = 1; i < N; i++) {
        c = A[i][i] + A[i][i - 1] * alpha[i - 1];
        alpha[i] = -A[i][i + 1] / c;
        beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / c;
    }

    x[N - 1] = beta[N - 1];

    for (int i = N - 2; i >= 0; i--) {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }
}


/**
 * Метод Холецкого
 */
void venediktovayap::lab4() {


    double **S = new double *[N];
    for (int i = 0; i < N; i++) {
        S[i] = new double[N];
        for (int j = 0; j < N; j++) S[i][j] = 0;
    }

    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if (j == i) {
                if (i == 0) {
                    S[i][i] = sqrt(A[0][0]);
                } else {
                    for (int k = 0; k < i; k++) A[i][i] -= S[k][i] * S[k][i];
                    S[i][i] = sqrt(A[i][i]);
                }
            } else {
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
void venediktovayap::lab5() {
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


double *MulVecToMatrix(int N, double *A[], double b[]) {
    double *temp = new double[N];
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
void venediktovayap::lab6() {

    double *r = new double[N];
    double *xk = new double[N];
    double eps = 1.e-17;
    double *x1 = b;
    double t;
    double maxDelta;
    do {
        double *temp = MulVecToMatrix(N, A, x1);

        for (int i = 0; i < N; i++) {
            r[i] = temp[i] - b[i];
        }
        double *Ar = MulVecToMatrix(N, A, r);

        double Scalar1, Scalar2;
        Scalar1 = ScalarMul(N, Ar, r);
        Scalar2 = ScalarMul(N, Ar, Ar);
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
void venediktovayap::lab7() {
}

double **TransposeMatrix(double **&m, int n) {
    double **temp = new double *[n];
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

double **MulMatrixes(double **&m1, double **&m2, int n) {
    double **temp = new double *[n];
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

void venediktovayap::lab8() {
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

        double **m_H = new double *[N];
        for (int i = 0; i < N; ++i) {
            m_H[i] = new double[N];
        }

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i == j) {
                    m_H[i][j] = 1;
                } else {
                    m_H[i][j] = 0;
                }
            }
        }

        m_H[maxi][maxj] = -sin(Phi);
        m_H[maxj][maxi] = sin(Phi);

        m_H[maxi][maxi] = cos(Phi);
        m_H[maxj][maxj] = cos(Phi);

        double **t_H = TransposeMatrix(m_H, N);

        double **t_A = MulMatrixes(t_H, A, N);

        A = MulMatrixes(t_A, m_H, N);
    } while (MaxEl >= eps);

    for (int i = 0; i < N; ++i) {
        cout << "lambda  = " << A[i][i] << endl;
    }
}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void venediktovayap::lab9() {
    double eps = 1.e-15;
    double lambda_k, lambda_k1 = 0;
    double *x_ = new double[N];
    double k;
    for (int i = 0; i < N; i++) {
        x[i] = b[i];
    }

    do {
        double sum1 = 0, sum2 = 0;
        lambda_k = lambda_k1;
        x_ = MulVecToMatrix(N, A, x);

        for (int i = 0; i < N; i++) {
            sum1 += x_[i];
            sum2 += x[i];
        }

        lambda_k1 = sum1 / sum2;
        k = sqrt(ScalarMul(N, x_, x_));
        for (int i = 0; i < N; i++) {
            x[i] = x_[i] / k;
        }

    } while (abs(lambda_k1 - lambda_k) >= eps);

    cout << "Max value of Lambda = " << lambda_k1 << endl;
    delete[] x_;
}


std::string venediktovayap::get_name() {
    return "R.V. Zhalnin";
}

#include "ryabikinks.h"

/**
 * Введение в дисциплину
 */
void ryabikinks::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void ryabikinks::lab2()
{
    int i = 0;
    while (i < N) {
        int maxColumn = i;
        for (int j = i + 1; j < N; j++) {
            if (abs(A[j][i]) > abs(A[maxColumn][i])) {
                maxColumn = j;
            }
        }
        if (maxColumn != i) {
            swap(A[i], A[maxColumn]);
            swap(b[i], b[maxColumn]);
        }
        for (int j = i + 1; j < N; j++) {
            double ratio = A[j][i] / A[i][i];
            for (int k = i; k < N; k++) {
                A[j][k] -= (ratio * A[i][k]);
            }
            b[j] -= (ratio * b[i]);
        }
        i++;
    }

    for (int i = N - 1; i >= 0; i--) {
        for (int j = i + 1; j < N; j++) {
            b[i] -= (A[i][j] * x[j]);
        }
        x[i] = b[i] / A[i][i];
    }
}


void fuller(double mass[], int N) {
	for (int i = 0; i < N; i++)
	{
		mass[i] = 0.0;
	}
}

/**
 * Метод прогонки
 */
void ryabikinks::lab3()
{
    double y = 0.0;
    double beta[N];
    double alpha[N];
    fuller(beta, N);
    fuller(alpha, N);

    for (int j = 0; j < N; j++)
    {
        if (j == 0) {
            y = A[j][j];
            alpha[j] = -A[j][j + 1] / y;
            beta[j] = b[j] / y;
            continue;
        }

        for (int i = 1; i < N - 1; i++)
        {
            y = A[i][i] + A[i][i - 1] * alpha[i - 1];
            alpha[i] = -A[i][i + 1] / y;
            beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / y;
        }

        if (j == N - 1) {
            y = A[j][j] + A[j][j - 1] * alpha[j - 1];
            beta[j] = (b[j] - A[j][j - 1] * beta[j - 1]) / y;
            continue;
        }
    }
    for (int i = N - 1; i >= 0; i--)
    {
        if (i == N - 1) {
            x[i] = beta[i];
            continue;
        }
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }

}

void MatrixTransposition(double** A, double** B, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            B[j][i] = A[i][j];
        }
    }
}
void Gauss(double** A, double* x, double* b, int N) {
    int i = 0;
    while (i < N) {
        int maxColumn = i;
        for (int j = i + 1; j < N; j++) {
            if (abs(A[j][i]) > abs(A[maxColumn][i])) {
                maxColumn = j;
            }
        }
        if (maxColumn != i) {
            swap(A[i], A[maxColumn]);
            swap(b[i], b[maxColumn]);
        }
        for (int j = i + 1; j < N; j++) {
            double ratio = A[j][i] / A[i][i];
            for (int k = i; k < N; k++) {
                A[j][k] -= (ratio * A[i][k]);
            }
            b[j] -= (ratio * b[i]);
        }
        i++;
    }

    for (int i = N - 1; i >= 0; i--) {
        for (int j = i + 1; j < N; j++) {
            b[i] -= (A[i][j] * x[j]);
        }
        x[i] = b[i] / A[i][i];
    }
}
void fuller(double** mass, int N) {
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            mass[i][j] = 0.0;
        }
    }
}


/**
 * Метод Холецкого
 */
void ryabikinks::lab4()
{
	double y[N];
    fuller(y, N);

    double** L = new double* [N];
    for (int i = 0; i < N; i++) {
        L[i] = new double[N];
    }
    fuller(L, N);

    double** Lt = new double* [N];
    for (int i = 0; i < N; i++) {
        Lt[i] = new double[N];
    }
    fuller(Lt, N);

    L[0][0] = sqrt(A[0][0]);
    for (int i = 1; i < N; i++)
    {
        L[i][0] = A[i][0] / L[0][0];
    }

    for (int i = 1; i < N; i++)
    {
        double temp = 0;
        for (int p = 0; p < i; p++)
        {
            temp += L[i][p] * L[i][p];
        }

        L[i][i] = sqrt(A[i][i] - temp);
        temp = 0;
        for (int j = i + 1; j < N; j++)
        {
            for (int p = 0; p < i; p++)
            {
                temp += L[i][p] * L[j][p];
            }
            L[j][i] = (A[j][i] - temp) / L[i][i];
        }
    }
    MatrixTransposition(L, Lt, N);
    Gauss(L, y, b, N);
    Gauss(Lt, x, y, N);
}



/**
 * Метод Якоби или Зейделя
 */
void ryabikinks::lab5()
{
    double eps = 1e-10;
    for (int i = 0; i < N; i++) {
        x[i] = 1;
    }
    double* buffer = new double[N];
    bool Convergence = true;
    while (Convergence) {
        Convergence = false;
        for (int i = 0; i < N; i++) {
            buffer[i] = b[i];
            for (int j = 0; j < N; j++) {
                if (j != i) {
                    buffer[i] -= A[i][j] * x[j];
                }
            }
            buffer[i] /= A[i][i];
            if (fabs(buffer[i] - x[i]) > eps) {
                Convergence = true;
            }
            x[i] = buffer[i];
        }
    }
}



/**
 * Метод минимальных невязок
 */
void ryabikinks::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void ryabikinks::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void ryabikinks::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void ryabikinks::lab9()
{

}


std::string ryabikinks::get_name()
{
  return "K.S.Ryabikin";
}

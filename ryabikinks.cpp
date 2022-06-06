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
    delete[] buffer;
}

double* MultMatrixVect(double** A, double b[], int N) {
    double* blank = new double[N];
    for (int i = 0; i < N; i++)
    {
        double temp = 0;
        for (int j = 0; j < N; j++)
        {
            temp += A[i][j] * b[j];
        }
        blank[i] = temp;
    }
    return blank;
}
double ScalarProduct(double a[], double b[], int N) {
    double blank = 0;
    for (int i = 0; i < N; i++) {
        blank += a[i] * b[i];
    }
    return blank;
}

/**
 * Метод минимальных невязок
 */
void ryabikinks::lab6()
{
    double eps = 1e-5;
    double r[N];
    double x1[N];
    double t = 0.0;
    double norm;
    fuller(r,N);
    fuller(x1,N);
    for (int i = 0; i < N; i++)
    {
        x1[i] = b[i];
        x[i] = b[i];
    }
    double* buffer = new double[N];
    buffer = MultMatrixVect(A, x, N);
    for (int i = 0; i < N; i++)
    {
        r[i] = buffer[i] - b[i];
    }
    double* bufferAr = new double[N];
    double* bufferAx = new double[N];
    do {
        bufferAr = MultMatrixVect(A, r, N);
        t = ScalarProduct(bufferAr, r, N) / ScalarProduct(bufferAr, bufferAr, N);
        for (int i = 0; i < N; i++)
        {
            x1[i] = x[i] - t * r[i];
        }
        bufferAx = MultMatrixVect(A, x1, N);
        for (int i = 0; i < N; i++)
        {
            r[i] = bufferAx[i] - b[i];
        }
        norm = fabs(x1[0] - x[0]);
        x[0] = x1[0];
        for (int i = 1; i < N; i++)
        {
            if (fabs(x1[i] - x[i]) > norm) {
                norm = fabs(x1[i] - x[i]);

            }
            x[i] = x1[i];
        }
    } while (norm > eps);

    delete[] buffer;
    delete[] bufferAr;
    delete[] bufferAx;
}



/**
 * Метод сопряженных градиентов
 */
void ryabikinks::lab7()
{
    double eps = 1e-5;
    double x1[N];
    double r[N];
    double r1[N];
    double p[N];
    fuller(x1,N);
    fuller(r,N);
    fuller(r1,N);
    fuller(p,N);
    double alpha, beta = 0.0;
    double norm;

    for (int i = 0; i < N; i++)
    {
        x1[i] = b[i];
        x[i] = b[i];
    }
    double* buffer = new double[N];
    buffer = MultMatrixVect(A, x, N);
    for (int i = 0; i < N; i++)
    {
        r[i] = b[i] - buffer[i];
    }
    double* bufferAp = new double[N];
    double* bufferAr = new double[N];
    for (int i = 0; i < N; i++)
    {
        p[i] = r[i];
    }
    do {
        for (int i = 0; i < N; i++)
        {
            r1[i] = r[i];
        }
        bufferAr = MultMatrixVect(A, r, N);		
        bufferAp = MultMatrixVect(A, p, N);
        alpha = ScalarProduct(r, r, N) / ScalarProduct(p, bufferAp, N);
        for (int i = 0; i < N; i++)
        {
            x1[i] = x[i] + alpha * p[i];
            r1[i] = r[i] - alpha * bufferAp[i];
        }
        beta = ScalarProduct(r1, r1, N) / ScalarProduct(r, r, N);
        for (int i = 0; i < N; i++)
        {
            p[i] = r1[i] + beta * p[i];
        }
        norm = fabs(x1[0] - x[0]);
        x[0] = x1[0];
        r[0] = r1[0];
        for (int i = 1; i < N; i++)
        {
            if (fabs(x1[i] - x[i]) > norm) {
                norm = fabs(x1[i] - x[i]);

            }
            x[i] = x1[i];
            r[i] = r1[i];
        }
    } while (norm >= eps);

    delete[] buffer;
    delete[] bufferAp;
    delete[] bufferAr;
}

void MultMatrix(double** A, double** B, double** C, double eps, int N) {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			C[i][j] = 0;
			for (int k = 0; k < N; ++k) {
				C[i][j] += A[i][k] * B[k][j];
				if (abs(C[i][j]) < eps) {
					C[i][j] = 0.0;
				}
			}
		}
	}
}
void BinaryMatrix(double** A, int N) {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			if (i == j) {
				A[i][j] = 1;
			}
			else {
				A[i][j] = 0;
			}
		}
	}
}
/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void ryabikinks::lab8()
{
	double eps = 1e-10;
	double epsilon = 1e-3;

	double maxEl;
	do {
		maxEl = 0.0;
		int i0 = 0, j0 = 0;
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				if (i == j) continue;
				if (abs(A[i][j]) > maxEl) {
					maxEl = abs(A[i][j]);
					i0 = i;
					j0 = j;
				}
			}
		}
		double phi = 0.5 * atan((2 * A[i0][j0]) / (A[i0][i0] - A[j0][j0]));
		double sinus = sin(phi);
		double cosine = cos(phi);
		double** buffer = new double*[N];
		for (int i = 0; i < N; ++i) {
			buffer[i] = new double[N];
		}

		BinaryMatrix(buffer, N);

		buffer[i0][i0] = cosine;
		buffer[i0][j0] = -sinus;
		buffer[j0][i0] = sinus;
		buffer[j0][j0] = cosine;

		double** bufferT = new double*[N];
		for (int i = 0; i < N; i++) {
			bufferT[i] = new double[N];
		}

		fuller(bufferT, N);
		MatrixTransposition(buffer, bufferT, N);

		double** IntermediateMatrix = new double*[N];
		for (int i = 0; i < N; i++) {
			IntermediateMatrix[i] = new double[N];
		}
		fuller(IntermediateMatrix, N);

		MultMatrix(bufferT, A, IntermediateMatrix, eps, N);
	
		MultMatrix(IntermediateMatrix, buffer, A, eps, N);

	} while (maxEl >= epsilon);
	for (int i = 0; i < N; ++i) {
		x[i] = A[i][i];
	}
}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void ryabikinks::lab9()
{
    double eps = 1e-17;
    double Lambda;
    double LambdaMax = 0;
    double* x1 = new double[N];
    double koef;
    for (int i = 0; i < N; i++)
    {
        x[i] = b[i];
    }

    do
    {
        Lambda = LambdaMax;
        x1 = MultMatrixVect(A, x, N);
        double sum1 = 0, sum2 = 0;
        for (int i = 0; i < N; i++)
        {
            sum1 += x1[i];
            sum2 += x[i];
        }

        LambdaMax = sum1 / sum2;
        koef = sqrt(ScalarProduct(x1, x1, N));
        for (int i = 0; i < N; i++)
        {
            x[i] = x1[i] / koef;
        }

    } while (abs(LambdaMax - Lambda) >= eps);

    delete[] x1;

    cout << "Max Lambda " << LambdaMax << endl;
}


std::string ryabikinks::get_name()
{
  return "K.S.Ryabikin";
}

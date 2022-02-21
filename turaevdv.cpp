#include "turaevdv.h"

using namespace std;

void methodOfGauss(double ** mas, int position, int size);
void swap(double ** arr, int size, int fromPosition, int toPosition);
double ** arr;

/**
 * Введение в дисциплину
 */
void turaevdv::lab1() {
    cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void turaevdv::lab2() {
    arr = createSLOUMatrix(N);
    methodOfGauss(arr, 0, N);
    for (int i = 0; i < N; ++i) {
        x[i] = arr[i][N];
    }
}

double **turaevdv::createSLOUMatrix(int size) {
    read_file();
    double **mas = new double *[size];
    for (int i = 0; i < size; ++i) {
        mas[i] = new double[size + 1];
    }
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            mas[i][j] = A[i][j];
        }
        mas[i][size] = b[i];
    }
    return mas;
}

void swapToLargest(double **mas, int position, int size) {
    int current = position;;
    int largestPosition = position;
    while (current < size) {
        if (mas[current][position] > mas[largestPosition][position] && mas[current][position] != 0) {
            largestPosition = current;
        }
        ++current;
    }
    if (largestPosition != position) {
        swap(mas, size, position, largestPosition);
    }
}

void swap(double **arr, int size, int fromPosition, int toPosition) {
    for (int i = fromPosition; i < size + 1; ++i) {
        arr[fromPosition][i] += arr[toPosition][i];
        arr[toPosition][i] = arr[toPosition][i] - arr[toPosition][i];
        arr[fromPosition][i] -= arr[toPosition][i];
    }
}

void methodOfGauss(double **mas, int position, int size) {

    swapToLargest(mas, position, size);
    for (int i = position + 1; i < size + 1; ++i) {
        mas[position][i] = mas[position][i] / mas[position][position];
    }
    mas[position][position] = 1;

    for (int i = position + 1; i < size; ++i) {
        for (int j = position + 1; j < size + 1; ++j) {
            mas[i][j] -= mas[position][j] * mas[i][position];
        }
        mas[i][position] = 0;
    }

    if (position + 1 < size) {
        methodOfGauss(mas, position + 1, size);
    }

    for (int i = position - 1; i >= 0; --i) {
        mas[i][size] -= mas[position][size] * mas[i][position];
        mas[i][position] = 0;
    }
}

/**
 * Метод прогонки
 */
void turaevdv::lab3() {
    b[0] /= A[0][0];
    A[0][0] = -A[0][1]/A[0][0];
    for (int i = 1; i < N - 1; ++i) {
        double y = A[i][i] + A[i][i-1]*A[i-1][i-1];
        A[i][i] = -A[i][i+1]/y;
        b[i] = (b[i] - A[i][i-1]*b[i-1])/y;
    }

    x[N-1] = (b[N-1] - A[N-1][N-2]*b[N-2]) / (A[N-1][N-1] + A[N-1][N-2] * A[N-2][N-2]);

    for (int i = N - 2; i >= 0; --i) {
        x[i] = x[i+1] * A[i][i] + b[i];
    }
}


void printMatrix(double **mas, int size);
void straightRunning(double **mas, int size);
void reverse(double **mas, int size);

/**
 * Метод Холецкого
 */
void turaevdv::lab4() {
    double **l = new double *[N];
    double **transposeL = new double *[N];
    for (int i = 0; i < N; ++i) {
        l[i] = new double[N+1];
        transposeL[i] = new double[N+1];
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N + 1; ++j) {
            l[i][j] = 0;
            transposeL[i][j] = 0;
        }
    }
    l[0][0] = sqrt(A[0][0]);
    transposeL[0][0] = l[0][0];
    for (int i = 1; i < N; ++i) {
        l[i][0] = A[i][0] / l[0][0];
        transposeL[0][i] = l[i][0];
    }
    for (int i = 1; i < N; ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += l[i][j] * l[i][j];
        }
        l[i][i] = sqrt(A[i][i] - sum);
        transposeL[i][i] = l[i][i];
        for (int j = i; j < N; ++j) {
            double sum = 0;
            for (int k = 0; k < i; ++k) {
                sum += l[i][k] * l[j][k];
            }
            l[j][i] = (A[j][i] - sum) / l[i][i];
            transposeL[i][j] = l[j][i];
        }
    }
    for (int i = 0; i < N; ++i) {
        l[i][N] = b[i];
    }
    straightRunning(l, N);
    for (int i = 0; i < N; ++i) {
        transposeL[i][N] = l[i][N];
    }
    reverse(transposeL, N);
    for (int i = 0; i < N; ++i) {
        x[i] = transposeL[i][N];
    }
}

void printMatrix(double **mas, int size) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            cout << mas[i][j] << " ";
        }
        cout << "\n";
    }
}

void straightRunning(double **mas, int size) {
    for (int i = 0; i < size; ++i) {
        mas[i][size] /= mas[i][i];
        for (int j = i + 1; j < size; ++j) {
            mas[j][size] = mas[j][size] - mas[j][i] * mas[i][size];
        }
    }
}

void reverse(double **mas, int size) {
    for (int i = size-1; i >= 0; --i) {
        mas[i][size] /= mas[i][i];
        for (int j = i - 1; j >= 0; --j) {
            mas[j][size] -= mas[i][size] * mas[j][i];
        }
    }
}

/**
 * Метод Якоби или Зейделя
 */
void turaevdv::lab5() {
    double norm = 0;
    double e = 1.e-20;
    for (int i = 0; i < N; ++i) {
        x[i] = b[i];
    }

    do {
        norm = 0;
        for(int i = 0; i < N; ++i) {
            double sum = 0;
            for (int j = 0; j < N; ++j) {
                if (i != j) {
                    sum += A[i][j] * x[j];
                }
            }
            sum = (b[i] - sum)/A[i][i];
            if (norm < (fabs(sum - x[i]))) {
                norm = (fabs(sum - x[i]));
            }
            x[i] = sum;
        }
    } while (norm >= e);
}


double* mulMatrixOnVector(double** matrix, double* vector, int size);
/**
 * Метод минимальных невязок
 */
void turaevdv::lab6() {
    double * nextVectorX = new double[N];

    for (int i = 0; i < N; ++i) {
        nextVectorX[i] = b[i];
    }

    double eps = 1.e-19;
    double *r = new double[N];
    double norm;

    do {
        double * mulMatrixAOnX = mulMatrixOnVector(A, x, N);

        for (int i = 0; i < N; ++i) {
            r[i] = mulMatrixAOnX[i] - b[i];
        }

        double * mulMatrixOnR = mulMatrixOnVector(A, r, N);

        double t = 0;
        double sum = 0;
        for (int i = 0; i < N; ++i) {
            t += mulMatrixOnR[i] * r[i];
            sum += pow(mulMatrixOnR[i], 2);
        }

        t /= sum;

        norm = 0;
        for (int i = 0; i < N; ++i) {
            x[i] = nextVectorX[i] - t * r[i];
            if (norm < fabs(nextVectorX[i] - x[i])) {
                norm = fabs(nextVectorX[i] - x[i]);
            }
            nextVectorX[i] = x[i];
        }
    } while (norm >= eps);
}

double* mulMatrixOnVector(double** matrix, double* vector, int size) {
    double* result = new double [size];
    for (int i = 0; i < size; ++i) {
        result[i] = 0;
        for (int j = 0; j < size; ++j) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
    return result;
}

double scalarProduct(double* v1, double* v2, int size) {
    double sum = 0;
    for (int i = 0; i < size; ++i) {
        sum += v1[i] * v2[i];
    }
    return sum;
}


/**
 * Метод сопряженных градиентов
 */
void turaevdv::lab7() {
    double * nextVectorX = new double[N];

    for (int i = 0; i < N; ++i) {
        x[i] = b[i];
    }

    double eps = 1.e-25;
    double * r = new double[N];
    double * z = new double[N];
    double norm;

    double * mulMatrixAOnX = mulMatrixOnVector(A, x, N);

    for (int i = 0; i < N; ++i) {
        r[i] = b[i] - mulMatrixAOnX[i];
        z[i] = r[i];
    }

    do {
        double * mulMatrixAOnZ = mulMatrixOnVector(A, z, N);

        double alpha = scalarProduct(r, r, N) / scalarProduct(mulMatrixAOnZ, z, N);

        double * nextVectorR = new double[N];

        for (int i = 0; i < N; ++i) {
            nextVectorX[i] = x[i] + alpha * z[i];
            nextVectorR[i] = r[i] - alpha * mulMatrixAOnZ[i];
        }

        double beta = scalarProduct(nextVectorR, nextVectorR, N) / scalarProduct(r, r, N);

        norm = 0;
        for (int i = 0; i < N; ++i) {
            z[i] = nextVectorR[i] + beta * z[i];
            r[i] = nextVectorR[i];
            if (norm < fabs(nextVectorX[i] - x[i])) {
                norm = fabs(nextVectorX[i] - x[i]);
            }
            x[i] = nextVectorX[i];
        }
    } while (norm >= eps);
}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void turaevdv::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void turaevdv::lab9()
{

}


std::string turaevdv::get_name()
{
  return "D.V. Turaev";
}

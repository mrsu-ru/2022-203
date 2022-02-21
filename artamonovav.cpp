#include "artamonovav.h"
#include <iostream>

/**
 * Введение в дисциплину
 */

using namespace std;

void artamonovav::lab1()
{
  cout << "Hello World!" << endl;
}


void refactorRow(double **&matrix, double *&parameters, int size, int numRow, double num) {
    for (int i = 0; i < size; ++i) {
        matrix[numRow][i] /= num;
    }

    parameters[numRow] /= num;
}

void subtract(double **&matrix, double *&parameters, int size, int numMainRow, int numSubtractedRow) {
    for (int i = 0; i < size; ++i) {
        matrix[numSubtractedRow][i] -= matrix[numMainRow][i];
    }

    parameters[numSubtractedRow] -= parameters[numMainRow];
}

void swapRows(double **&matrix, double *&parameters, int size, int numRow1, int numRow2) {
    for (int i = 0; i < size; ++i) {
        double d = matrix[numRow1][i];
        matrix[numRow1][i] = matrix[numRow2][i];
        matrix[numRow2][i] = d;
    }

    double d = parameters[numRow1];
    parameters[numRow1] = parameters[numRow2];
    parameters[numRow2] = d;
}

int searchRow(double **&matrix, int size, int numBeginRow, int numColumn) {
    for (int i = numBeginRow; i < size; ++i) {
        if (matrix[i][numColumn] != 0) {
            return i;
        }
    }

    return -1;
}

void firstPartOfGaussAlgorithm(double **&matrix, double *&parameters, int size, int pos) {
    int t = pos;
    if (matrix[t][t] == 0) {
        t = searchRow(matrix, size, t, pos);

        if (t != -1) {
            swapRows(matrix, parameters, size, t, pos);
        }
    }

    if ((matrix[pos][pos] != 1) && (matrix[pos][pos] != 0)) {
        refactorRow(matrix, parameters, size, pos, matrix[pos][pos]);
    }

    if (pos != size - 1) {
        for (int i = pos + 1; i < size; ++i) {
            double d = matrix[i][pos];

            if (d != 0) {
                refactorRow(matrix, parameters, size, i, d);
                subtract(matrix, parameters, size, pos, i);
                refactorRow(matrix, parameters, size, i, 1 / d);
            }
        }
    }
}

void lastPartOfGaussAlgorithm(double **&matrix, double *&parameters, int size, int pos) {
    if (matrix[pos][pos] != 1) {
        refactorRow(matrix, parameters, size, pos, matrix[pos][pos]);
    }

    if (pos != 0) {
        for (int i = pos - 1; i >= 0; --i) {
            double d = matrix[i][pos];

            if (d != 0) {
                refactorRow(matrix, parameters, size, i, d);
                subtract(matrix, parameters, size, pos, i);
                refactorRow(matrix, parameters, size, i, 1 / d);
            }
        }
    }
}

/**
 * Метод Гаусса с выбором главного элемента
 */
void artamonovav::lab2() {
    for (int i = 0; i < N; ++i) {
        firstPartOfGaussAlgorithm(A, b, N, i);
    }

    for (int i = N - 1; i >= 0; --i) {
        lastPartOfGaussAlgorithm(A, b, N, i);
    }

    for (int i = 0; i < N; ++i) {
        x[i] = b[i];
    }
}

/**
 * Метод прогонки
 */
void artamonovav::lab3() {
    for (int i = 1; i < N; ++i) {
        b[i] -= A[i][i - 1] / A[i - 1][i - 1] * b[i - 1];
        A[i][i] -= A[i][i - 1] / A[i - 1][i - 1] * A[i - 1][i];
    }

    x[N - 1] = b[N - 1] / A[N - 1][N - 1];
    for (int i = N - 1; i >= 0; --i) {
        x[i] = (b[i] - A[i][i + 1] * x[i + 1]) / A[i][i];
    }
}

/**
 * Метод Холецкого
 */

void artamonovav::lab4() {
    double **matrixL = new double*[N];
    double **transposeMatrixL = new double*[N];

    for (int i = 0; i < N; ++i) {
        matrixL[i] = new double[N];
        transposeMatrixL[i] = new double[N];
    }

    matrixL[0][0] = sqrt(A[0][0]);
    transposeMatrixL[0][0] = matrixL[0][0];
    for (int i = 1; i < N; ++i) {
        matrixL[i][0] = A[i][0] / matrixL[0][0];
        transposeMatrixL[0][i] = matrixL[i][0];
    }

    for (int i = 1; i < N; ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += matrixL[i][j] * matrixL[i][j];
        }
        matrixL[i][i] = sqrt(A[i][i] - sum);
        transposeMatrixL[i][i] = matrixL[i][i];

        for (int j = i; j < N; ++j) {
            sum = 0;
            for (int k = 0; k < i; ++k) {
                sum += matrixL[i][k] * matrixL[j][k];
            }
            matrixL[j][i] = (A[j][i] - sum) / matrixL[i][i];
            transposeMatrixL[i][j] = matrixL[j][i];
        }
    }

    for (int i = 0; i < N; ++i) {
        firstPartOfGaussAlgorithm(matrixL, b, N, i);
    }

    for (int i = N - 1; i >= 0; --i) {
        lastPartOfGaussAlgorithm(transposeMatrixL, b, N, i);
    }

    for (int i = 0; i < N; ++i) {
        x[i] = b[i];
    }
}

/**
 * Метод Якоби или Зейделя
 */
void artamonovav::lab5() {
    double * x1 = new double[N];

    for (int i = 0; i < N; ++i) {
        x1[i] = b[i];
    }

    double eps = 1.e-25;
    double norm = 0;

    do {
        for (int i = 0; i < N; ++i) {
            x[i] = x1[i];
        }

        for (int i = 0; i < N; ++i) {
            double sum = 0;
            for (int j = 0; j < N; ++j) {
                if (i != j) {
                    sum += A[i][j] * x1[j];
                }
            }
            x[i] = (b[i] - sum) / A[i][i];
        }

        norm = fabs(x1[0] - x[0]);
        for (int i = 1; i < N; ++i) {
            if (norm < fabs(x1[i] - x[i])) {
                norm = fabs(x1[i] - x[i]);
            }
        }

        for (int i = 0; i < N; ++i) {
            x1[i] = x[i];
        }
    } while (norm >= eps);
}


double * getMulMatrixOnVector(double **&matrix, double *&vector, int size) {
    double * result = new double[size];

    for (int i = 0; i < size; ++i) {
        result[i] = 0;
        for (int j = 0; j < size; ++j) {
            result[i] += matrix[i][j] * vector[j];
        }
    }

    return result;
}

/**
 * Метод минимальных невязок
 */
void artamonovav::lab6() {
    double * x1 = new double[N];

    for (int i = 0; i < N; ++i) {
        x1[i] = b[i];
    }

    double eps = 1.e-19;
    double *r = new double[N];
    double norm;

    do {
        for (int i = 0; i < N; ++i) {
            x[i] = x1[i];
        }

        double * A_x = getMulMatrixOnVector(A, x, N);

        for (int i = 0; i < N; ++i) {
            r[i] = A_x[i] - b[i];
        }

        double * A_r = getMulMatrixOnVector(A, r, N);

        double t = 0;
        double sum = 0;
        for (int i = 0; i < N; ++i) {
            t += A_r[i] * r[i];
            sum += pow(A_r[i], 2);
        }

        t /= sum;

        for (int i = 0; i < N; ++i) {
            x[i] = x1[i] - t * r[i];
        }

        norm = fabs(x1[0] - x[0]);
        for (int i = 1; i < N; ++i) {
            if (norm < fabs(x1[i] - x[i])) {
                norm = fabs(x1[i] - x[i]);
            }
        }

        for (int i = 0; i < N; ++i) {
            x1[i] = x[i];
        }
    } while (norm >= eps);
}


double getScalarMul(double *&vec1, double *&vec2, int count) {
    double result = 0;

    for (int i = 0; i < count; ++i) {
        result += vec1[i] * vec2[i];
    }

    return result;
}

/**
 * Метод сопряженных градиентов
 */
void artamonovav::lab7() {
    double * x1 = new double[N];

    for (int i = 0; i < N; ++i) {
        x1[i] = b[i];
    }

    double eps = 1.e-25;
    double * r = new double[N];
    double * z = new double[N];
    double norm;

    double * A_x = getMulMatrixOnVector(A, x1, N);

    for (int i = 0; i < N; ++i) {
        r[i] = b[i] - A_x[i];
        z[i] = r[i];
    }

    do {
        double * A_z = getMulMatrixOnVector(A, z, N);

        double alpha = getScalarMul(r, r, N) / getScalarMul(A_z, z, N);

        double * r1 = new double[N];

        for (int i = 0; i < N; ++i) {
            x[i] = x1[i] + alpha * z[i];
            r1[i] = r[i] - alpha * A_z[i];
        }

        double beta = getScalarMul(r1, r1, N) / getScalarMul(r, r, N);

        for (int i = 0; i < N; ++i) {
            z[i] = r1[i] + beta * z[i];
            r[i] = r1[i];
        }

        norm = fabs(x1[0] - x[0]);
        for (int i = 1; i < N; ++i) {
            if (norm < fabs(x1[i] - x[i])) {
                norm = fabs(x1[i] - x[i]);
            }
        }

        for (int i = 0; i < N; ++i) {
            x1[i] = x[i];
        }
    } while (norm >= eps);
}


double ** getMulOfMatrixes(double **&matrix1, double **&matrix2, int count) {
    double ** result = new double * [count];
    for (int i = 0; i < count; ++i) {
        result[i] = new double[count];
    }

    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < count; ++j) {
            result[i][j] = 0;
            for (int k = 0; k < count; ++k) {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }

    return result;
}

double ** getTransposeMatrix(double **&matrix, int count) {
    double ** result = new double * [count];
    for (int i = 0; i < count; ++i) {
        result[i] = new double[count];
    }

    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < count; ++j) {
            result[i][j] = matrix[j][i];
        }
    }

    return result;
}

void printMatrix(double **&matrix, int count) {
    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < count; ++j) {
            cout << matrix[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;
    cout << endl;
}

/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void artamonovav::lab8() {
    double eps = 1.e-1;
    double norm;

    do {
        norm = 0;
        int iOfMaxElem = 0, jOfMaxElem = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i < j) {
                    if (norm < fabs(A[i][j])) {
                        norm = fabs(A[i][j]);
                        iOfMaxElem = i;
                        jOfMaxElem = j;
                    }
                }
            }
        }

        double sinPhi, cosPhi;

        if (A[iOfMaxElem][iOfMaxElem] != A[jOfMaxElem][jOfMaxElem]) {
            double p = 2 * A[iOfMaxElem][jOfMaxElem] / (A[iOfMaxElem][iOfMaxElem] - A[jOfMaxElem][jOfMaxElem]);

            double p1 = pow(1 + p * p, 0.5);

            sinPhi = ((p < 0) ? -1 : 1) * pow(0.5 * (1 - 1 / p1), 0.5);
            cosPhi = pow(0.5 * (1 + 1 / p1), 0.5);
        } else {
            sinPhi = sin(M_PI / 4);
            cosPhi = cos(M_PI / 4);
        }

        double ** matrixH = new double * [N];
        for (int i = 0; i < N; ++i) {
            matrixH[i] = new double[N];
        }

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i == j) {
                    matrixH[i][j] = 1;
                } else {
                    matrixH[i][j] = 0;
                }
            }
        }

        matrixH[iOfMaxElem][jOfMaxElem] = -sinPhi;
        matrixH[jOfMaxElem][iOfMaxElem] = sinPhi;

        matrixH[iOfMaxElem][iOfMaxElem] = cosPhi;
        matrixH[jOfMaxElem][jOfMaxElem] = cosPhi;

        double ** transposeH = getTransposeMatrix(matrixH, N);

        double ** transposeHMultipliedByA = getMulOfMatrixes(transposeH, A, N);

        A = getMulOfMatrixes(transposeHMultipliedByA, matrixH, N);
    } while (norm >= eps);

    for (int i = 0; i < N; ++i) {
        cout << A[i][i] << endl;
    }
}

/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void artamonovav::lab9() {
    double maxEigenvalueModulo;
    double num = 0;
    double eps = 1.e-10;

    double * x1 = new double[N];
    for (int i = 0; i < N; ++i) {
        x1[i] = b[i];
    }

    do {
        maxEigenvalueModulo = num;

        for (int i = 0; i < N; ++i) {
            x[i] = x1[i];
        }

        x1 = getMulMatrixOnVector(A, x, N);

        double sum1 = 0, sum2 = 0;
        for (int i = 0; i < N; ++i) {
            sum1 += x1[i];
            sum2 += x[i];
        }
        num = sum1 / sum2;

        double r = pow(getScalarMul(x1, x1, N), 0.5);
        for (int i = 0; i < N; ++i) {
            x1[i] /= r;
        }
    } while (fabs(num - maxEigenvalueModulo) >= eps);

    cout << "Max Eigenvalue Modulo: " << maxEigenvalueModulo << endl;
}

std::string artamonovav::get_name() {
  return "Artamonov Alexey";
}

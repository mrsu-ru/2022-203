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

void firstPart(double **&matrix, double *&parameters, int size, int pos) {
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

void lastPart(double **&matrix, double *&parameters, int size, int pos) {
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
void artamonovav::lab2()
{
    for (int i = 0; i < N; ++i) {
        firstPart(A, b, N, i);
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << A[i][j] << "\t";
        }
        cout << endl;
    }

    for (int i = N - 1; i >= 0; --i) {
        lastPart(A, b, N, i);
    }

    for (int i = 0; i < N; ++i) {
        x[i] = b[i];
    }

}



/**
 * Метод прогонки
 */
void artamonovav::lab3()
{
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
void artamonovav::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void artamonovav::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void artamonovav::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void artamonovav::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void artamonovav::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void artamonovav::lab9()
{

}


std::string artamonovav::get_name()
{
  return "Artamonov Alexey";
}

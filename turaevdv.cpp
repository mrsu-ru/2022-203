#include "turaevdv.h"

using namespace std;

void methodOfGauss(double **mas, int position, int size);

void swap(double **arr, int size, int fromPosition, int toPosition);

double **arr;

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

}


/**
 * Метод минимальных невязок
 */
void turaevdv::lab6() {

}


/**
 * Метод сопряженных градиентов
 */
void turaevdv::lab7() {

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void turaevdv::lab8() {

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void turaevdv::lab9() {

}


std::string turaevdv::get_name() {
    return "D.V. Turaev";
}

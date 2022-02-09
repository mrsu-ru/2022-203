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

double ** turaevdv::createSLOUMatrix(int size) {
    read_file();
    double ** mas = new double* [size];
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

void swapToLargest(double ** mas, int position, int size) {
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

void swap(double ** arr, int size, int fromPosition, int toPosition) {
    for (int i = fromPosition; i < size + 1; ++i) {
        arr[fromPosition][i] += arr[toPosition][i];
        arr[toPosition][i] = arr[toPosition][i] - arr[toPosition][i];
        arr[fromPosition][i] -= arr[toPosition][i];
    }
}

void methodOfGauss(double ** mas, int position, int size) {

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
void turaevdv::lab3()
{

}



/**
 * Метод Холецкого
 */
void turaevdv::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void turaevdv::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void turaevdv::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void turaevdv::lab7()
{

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

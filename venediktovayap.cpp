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

}


/**
 * Метод Холецкого
 */
void venediktovayap::lab4() {

}


/**
 * Метод Якоби или Зейделя
 */
void venediktovayap::lab5() {

}


/**
 * Метод минимальных невязок
 */
void venediktovayap::lab6() {

}


/**
 * Метод сопряженных градиентов
 */
void venediktovayap::lab7() {

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void venediktovayap::lab8() {

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void venediktovayap::lab9() {

}


std::string venediktovayap::get_name() {
    return "R.V. Zhalnin";
}

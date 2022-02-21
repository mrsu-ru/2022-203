#include "fedinda.h"
#include <iostream>
/**
 * Введение в дисциплину
 */
using namespace std;

void fedinda::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void fedinda::lab2()
{
  for (int i = 0; i < N; i++) {
		int max = i;
		for (int j = i; j < N; j++) {
			if (A[j][i] > A[max][i]) {
				max = j;
			}
		}
		if (max != i) {
			swap(A[i], A[max]);
			swap(b[i], b[max]);
		}
		double mainElem = A[i][i];
		for (int j = i; j < N; j++) {
			A[i][j] /= mainElem;
		}
		b[i] /= mainElem;

		for (int j = 0; j < N; j++) {
			if (j != i) {
				mainElem = A[j][i];
				for (int k = 0; k < N;k++) {
					A[j][k] -= mainElem * A[i][k];
				}
				b[j] -= mainElem * b[i];
			}
		}
	}
  
	for (int i = 0; i < N;i++) {
		x[i] = b[i];
	}
}



/**
 * Метод прогонки
 */
void fedinda::lab3()
{

}



/**
 * Метод Холецкого
 */
void fedinda::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void fedinda::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void fedinda::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void fedinda::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void fedinda::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void fedinda::lab9()
{

}


std::string fedinda::get_name()
{
  return "D.A. Fedin";
}

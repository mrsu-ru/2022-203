#include "akimovada.h"

/**
 * Введение в дисциплину
 */
void akimovada::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void akimovada::lab2()
{
    int position = 0;
	for (int i = 0; i < N; i++)
	{
		int k = i;
		for (int q = i + 1; q < N; q++)
		{
			if (A[k][i] < A[q][i] && A[k][i] == 0)
			{
				k = q;
			}
		}

		if (k != i){
			for (int j = 0; j < N; j++)
			{
				swap(A[i][j], A[k][j]);
			}
		}

		double main_element = A[i][i];
		for (int j = 0; j < N; j++)
		{
			A[i][j] /= main_element;
		}

		b[i] /= main_element;
		for (int j = 0; j < N; j++)
		{
			if (j != i)
			{
				double main_div = -A[j][i];
				for (int z = 0; z < N; z++)
				{
					A[j][z] += main_div * A[i][z];
				}

				b[j] += main_div * b[i];
			}
		}
	}

	for (int i = 0; i < N; i++) {
		x[i] = b[i];
	}

}



/**
 * Метод прогонки
 */
void akimovada::lab3()
{

}



/**
 * Метод Холецкого
 */
void akimovada::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void akimovada::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void akimovada::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void akimovada::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void akimovada::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void akimovada::lab9()
{

}


std::string akimovada::get_name()
{
  return "D. A. Akimova";
}

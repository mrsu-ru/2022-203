#include "nuyanzinma.h"

/**
 * Введение в дисциплину
 */
void nuyanzinma::lab1()
{
  cout << "Hello World!" << endl;
}
void printMatrix(int N, double**& A)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
}

void swapRows(int i, int j, int N, double**& A, double* &b)
{
	for (int k = 0; k < N; k++)
	{
		double c = A[i][k];
		A[i][k] = A[j][k];
		A[j][k] = c;
	}

	double c = b[i];
	b[i] = b[j];
	b[j] = c;
}

int findRowWithMaxFirstNumber(int i, int j, int N, double**& A)
{
	int n = i;
	for (int k = i; k < N; k++)
	{
		if (A[k][j] > A[n][j])
		{
			n = k;
		}
	}
	return n;
}

/**
 * Метод Гаусса с выбором главного элемента
 */
void nuyanzinma::lab2()
{
	// Прямой ход метода
	for (int i = 0; i < N; i++)
	{
		int rowWithMaxNumber = findRowWithMaxFirstNumber(i, i, N, A);
		swapRows(i, rowWithMaxNumber, N, A, b);
		double firstNumber = A[i][i];
		for (int k = i; k < N; k++) // Деление текущей строки на первый элемент этой строки
		{
			A[i][k] /= firstNumber;
		}
		b[i] /= firstNumber;

		for (int j = i + 1; j < N; j++)
		{
			double multiplier = A[j][i];
			for (int k = i; k < N; k++)
			{
				A[j][k] -= A[i][k] * multiplier;
			}
			b[j] -= b[i] * multiplier;
		}
	}
	
	// Обратный ход метода
	for (int i = N - 1; i > 0; i--)
	{
		for (int j = i - 1; j >= 0; j--)
		{
			double multiplier = A[j][i];
			for (int k = i; k >= 0; k--)
			{
				A[j][k] -= A[i][k] * multiplier;
			}
			b[j] -= b[i] * multiplier;
		}
	}
	for (int i = 0; i < N; i++)
	{
		x[i] = b[i];
	}
}



/**
 * Метод прогонки
 */
void nuyanzinma::lab3()
{
	double* Vk = new double[N];
	double* Uk = new double[N];

	Vk[0] = A[0][1] / -A[0][0];
	Uk[0] = b[0] / A[0][0];

	for (int i = 1; i < N - 1; i++)
	{
		Vk[i] = A[i][i + 1] / (-A[i][i] - A[i][i - 1] * Vk[i - 1]);
		Uk[i] = (A[i][i - 1] * Uk[i - 1] - b[i]) / (-A[i][i] - A[i][i - 1] * Vk[i - 1]);
	}

	Vk[N - 1] = 0;
	Uk[N - 1] = (A[N - 1][N - 2] * Uk[N - 2] - b[N - 1]) / (-A[N - 1][N - 1] - A[N - 1][N - 2] * Vk[N - 2]);

	x[N - 1] = Uk[N - 1];
	for (int i = N - 1; i > 0; i--)
	{
		x[i - 1] = Vk[i - 1] * x[i] + Uk[i - 1];
	}
}



/**
 * Метод Холецкого
 */
void nuyanzinma::lab4()
{

}



/**
 * Метод Якоби или Зейделя
 */
void nuyanzinma::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void nuyanzinma::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void nuyanzinma::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void nuyanzinma::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void nuyanzinma::lab9()
{

}


std::string nuyanzinma::get_name()
{
  return "M.A. Nuyanzin";
}
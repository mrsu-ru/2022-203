#include "negryame.h"

/**
 * Введение в дисциплину
 */
void negryame::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void negryame::lab2()
{
	for (int i = 0; i < N; i++) {

		int maxLine = i;

		for (int j = i + 1; j < N; j++) {
			if (abs(A[j][i]) > abs(A[maxLine][i])) {
				maxLine = j;
			}
		}

		if (maxLine != i) {
			for (int j = 0; j < N; j++) {
				swap(A[maxLine][j], A[i][j]);
				swap(b[maxLine], b[i]);
			}
		}

		for (int j = i + 1; j < N; j++) {
			double k = A[j][i] / A[i][i];
			for (int c = i; c < N; c++) {
				A[j][c] -= k * A[i][c];
			}
			b[j] -= (k * b[i]);
		}
	}
	for (int i = N - 1; i >= 0; i--) {
		for (int j = i + 1; j < N; j++) {
			b[i] -= A[i][j] * x[j];
		}
		x[i] = b[i] / A[i][i];
	}

}



/**
 * Метод прогонки
 */
void negryame::lab3()
{
	double* P = new double[N]; 
	double* Q = new double[N]; 

	P[0] = -A[0][1] / A[0][0];
	Q[0] = b[0] / A[0][0];

	for (int i = 1; i < N; i++) 
	{
		P[i] = A[i][i + 1] / (-A[i][i] - A[i][i - 1] * P[i - 1]);
		Q[i] = (-b[i] + A[i][i - 1] * Q[i - 1]) / (-A[i][i] - A[i][i - 1] * P[i - 1]);
	}

	x[N - 1] = Q[N - 1];
	for (int i = N - 2; i >= 0; i--) 
		x[i] = P[i] * x[i + 1] + Q[i];

	delete[] P;
	delete[] Q;
	

}





/**
 * Метод Холецкого
 */
void negryame::lab4()
{
	double** S = new double* [N];
	for (int i = 0; i < N; i++) {
		S[i] = new double[N];
		for (int j = 0; j < N; j++)
			S[i][j] = 0.0;
	}
	double* d = new double[N];


	for (int i = 0; i < N; i++)
	{
		double temp = A[i][i];
		for (int e = 0; e <= i - 1; e++) {
			temp -= d[e] * S[e][i] * S[e][i];
		}
		d[i] = signbit(temp) == false ? 1 : -1;
		S[i][i] = sqrt(temp * d[i]);


		for (int j = i + 1; j < N; j++)
		{
			double t = 0;
			for (int k = 0; k <= j - 1; k++)
				t += d[k] * S[k][i] * S[k][j];
			S[i][j] = (A[i][j] - t) / (d[i] * S[i][i]);
		}
	}


	for (int i = 0; i < N; i++) {
		b[i] /= S[i][i];
		for (int j = i + 1; j < N; j++)
			b[j] -= b[i] * S[i][j];
	}

	for (int i = 0; i < N; i++)
		for (int j = i; j < N; j++)
			S[i][j] *= d[i];

	for (int i = N - 1; i >= 0; i--) {
		b[i] /= S[i][i];
		for (int j = i - 1; j >= 0; j--)
			b[j] -= b[i] * S[j][i];
	}

	x = b;
}


}



/**
 * Метод Якоби или Зейделя
 */
void negryame::lab5()
{

	double* new_x = new double[N],
		eps = 1.e-10;
	bool condition;

	for (int i = 0; i < N; i++)
		x[i] = 1;

	do
	{
		condition = false;
		for (int i = 0; i < N; i++)
		{
			new_x[i] = b[i];
			for (int j = 0; j < N; j++)
			{
				if (i == j) continue;
				new_x[i] -= A[i][j] * x[j];
			}

			new_x[i] /= A[i][i];
			if (!condition) condition = (fabs(new_x[i] - x[i]) > eps);
			x[i] = new_x[i];
		}

	} while (condition);

	delete[] new_x;
}



/**
 * Метод минимальных невязок
 */
void negryame::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void negryame::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void negryame::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void negryame::lab9()
{

}


std::string negryame::get_name()
{
  return "M.E. Negrya";
}

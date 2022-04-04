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
	for (int i = 0; i < N; i++)
	{
		int k = i;
		for (int q = i + 1; q < N; q++)
		{
			if (abs(A[k][i]) < abs(A[q][i]) || A[k][i] == 0)
			{
				k = q;
			}
		}

		if (k != i){
			for (int j = 0; j < N; j++)
			{
				swap(A[i][j], A[k][j]);
                swap(b[i], b[k]);
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

void Gauss(double **A, double *b, int N)
{
    for (int i = 0; i < N; i++)
    {
        int k = i;
        for (int q = i + 1; q < N; q++)
        {
            if (abs(A[k][i]) < abs(A[q][i]) || A[k][i] == 0)
            {
                k = q;
            }
        }

        if (k != i){
            for (int j = 0; j < N; j++)
            {
                swap(A[i][j], A[k][j]);
                swap(b[i], b[k]);
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
}

void printM(double**matrix, int N)
{
    for (int i = 0; i < N; i++)
    {
        cout<<endl;
        for (int j = 0; j < N; j++)
        {
            cout<<matrix[i][j]<<" ";
        }
    }

    cout<<endl;
}
void printA(double* array, int N)
{
    cout<<endl;
    for (int i = 0; i < N; i++)
    {
        cout<<array[i]<<" ";
    }
    cout<<endl;
}



/**
 * Метод прогонки
 */
void akimovada::lab3()
{
	double alpha[N - 1];
	double beta[N];

	double y = A[0][0];
	alpha[0] = -A[0][1] / y;
	beta[0] = b[0] / y;

	for (int i = 1; i < N - 1; i++)
	{
		y = A[i][i] + A[i][i - 1] * alpha[i - 1];
		alpha[i] = -A[i][i + 1] / y;
		beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / y;
	}

	beta[N - 1] = (b[N - 1] - A[N - 1][N - 2] * beta[N - 2]) / (A[N - 1][N - 1] + A[N - 1][N - 2] * alpha[N - 2]);
	x[N - 1] = beta[N - 1];

	for (int i = N - 2; i >= 0; i--)
	{
		x[i] = alpha[i] * x[i + 1] + beta[i];
	}

}



/**
 * Метод Холецкого
 */
void akimovada::lab4()
{
    double**L = new double *[N];
    for (int i = 0; i < N; i++)
    {
        L[i] = new double[N];
        for (int j = 0; j < N; j++)
        {
            L[i][j] = 0;
        }
    }

    L[0][0] = sqrt(A[0][0]);
    for (int i = 1; i < N; i++)
    {
        L[i][0] = A[i][0] / L[0][0];
    }

    for (int i = 1; i < N; i++)
    {
        double temp = 0;
        for (int p = 0; p < i; p++)
        {
            temp += L[i][p] * L[i][p];
        }

        L[i][i] = sqrt(A[i][i] - temp);
        temp = 0;
        for (int j = i + 1; j < N; j++)
        {
            for (int p = 1; p < i - 1; p++)
            {
                temp += L[i][p] * L[j][p];
            }
            L[j][i] = (A[j][i] - temp) / L[i][i];
        }
    }

    double **L_T = new double*[N];

    for (int i = 0; i < N; i++)
    {
        L_T[i] = new double[N];
        for (int j = 0; j < N; j++)
        {
            L_T[i][j] = L[j][i];
        }
    }

    Gauss(L, b, N);
    Gauss(L_T, b, N);

    for (int i = 0; i < N; i++)
    {
        x[i] = b[i];
    }


}



/**
 * Метод Якоби или Зейделя
 */
void akimovada::lab5()
{
    double **B = new double*[N];
    for (int i = 0; i < N; i++)
    {
        B[i] = new double[N];
        for (int j = 0; j < N + 1; j++)
        {
            if (i == j)
                B[i][j] = 0;
            else
                B[i][j] = -A[i][j] / A[i][i];
        }
    }
    double* d = new double[N];
    double* dconst = new double[N];
    for (int i = 0; i < N; i++)
    {
        d[i] = b[i] / A[i][i];
        dconst[i] = b[i] / A[i][i];
    }

    double norm;
    do {
        for (int i = 0; i < N; i++)
        {
            double temp = 0;
            for (int j = 0; j < N; j++)
            {
                temp += B[i][j] * d[j];
            }

            x[i] = temp + dconst[i];
        }
        norm = fabs(d[0] - x[0]);
        for (int i = 1; i < N; i++)
        {
            if (fabs(d[i] - x[i]) < norm)
            {
                norm = fabs(d[i] - x[i]);
            }

            d[i] = x[i];
        }

    } while (norm > 1.e-18);
}



/**
 * Метод минимальных невязок
 */
void akimovada::lab6()
{
    double eps = 1.e-18, *x_0 = new double[N], *r_0 = new double[N], *temp = new double[N], modul = 0;
    for (int i = 0; i < N; i++)
    {
        x_0[i] = b[i];
        temp[i] = 0;
        for (int j = 0; j < N; j++)
        {
            temp[i] += x_0[j]*A[i][j];
        }

        r_0[i] = b[i] - temp[i];
        modul += r_0[i] * r_0[i];
    }
    
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

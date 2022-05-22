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

void forwardGaussMethod(int N, double**& L, double*& y, double*& b)
{
	for (int i = 0; i < N; i++)
	{
		double firstNumber = L[i][i];
		for (int k = i; k < N; k++)
		{
			L[i][k] /= firstNumber;
		}
		b[i] /= firstNumber;

		for (int j = i + 1; j < N; j++)
		{
			double multiplier = L[j][i];
			for (int k = i; k < N; k++)
			{
				L[j][k] -= L[i][k] * multiplier;
			}
			b[j] -= b[i] * multiplier;
		}
	}
	for (int i = 0; i < N; i++)
	{
		y[i] = b[i];
	}
}

void reverseGaussMethod(int N, double**& Ltr, double*& x, double*& y)
{
	for (int i = N - 1; i >= 0; i--)
	{
		double ltrDiag = Ltr[i][i];
		for (int k = i; k < N; k++)
		{
			Ltr[i][k] /= ltrDiag;
		}
		y[i] /= ltrDiag;

		for (int j = i - 1; j >= 0; j--)
		{
			double multiplier = Ltr[j][i];
			for (int k = i; k >= 0; k--)
			{
				Ltr[j][k] -= Ltr[i][k] * multiplier;
			}
			y[j] -= y[i] * multiplier;
		}
	}
	for (int i = 0; i < N; i++)
	{
		x[i] = y[i];
	}
}

/**
 * Метод Холецкого
 */
void nuyanzinma::lab4()
{
	double** L = new double* [N];
	for (int i = 0; i < N; i++)
	{
		L[i] = new double[N];
	}
	for (int i = 0; i < N; i++)
	{
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
		double sum = 0;
		for (int p = 0; p < i; p++)
		{
			sum += L[i][p] * L[i][p];
		}
		L[i][i] = sqrt(A[i][i] - sum);

		for (int j = i + 1; j < N; j++)
		{
			double sum = 0;
			for (int p = 0; p < i; p++)
			{
				sum += L[i][p] * L[j][p];
			}
			L[j][i] = 1 / L[i][i] * (A[j][i] - sum);
		}
	}

	double** Ltr = new double* [N];
	for (int i = 0; i < N; i++)
	{
		Ltr[i] = new double[N];
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Ltr[i][j] = L[j][i];
		}
	}
	double* y = new double[N];
	forwardGaussMethod(N, L, y, b);
	reverseGaussMethod(N, Ltr, x, y);
}

double** productOfMatricies(int N, double**& A, double**& B)
{
	double** result = new double*[N];
	for (int i = 0; i < N; i++)
	{
		result[i] = new double[N];
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			double sum = 0.0;
			for (int k = 0; k < N; k++)
			{
				sum += A[i][k] * B[k][j];
			}
			result[i][j] = sum;
		}
	}

	return result;
}

double* productOfMatrixAndVector(int N, double**& A, double*& d)
{
	double* result = new double[N];
	for (int i = 0; i < N; i++)
	{
		double sum = 0;
		for (int j = 0; j < N; j++)
		{
			sum += A[i][j] * d[j];
		}
		result[i] = sum;
	}
	return result;
}

/**
 * Метод Якоби или Зейделя
 */
void nuyanzinma::lab5()
{
	//Метод Якоби
	double** D = new double* [N];
	double** Drev = new double* [N];
	double** DminusA = new double* [N];
	for (int i = 0; i < N; i++)
	{
		D[i] = new double[N];
		Drev[i] = new double[N];
		DminusA[i] = new double[N];
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			D[i][j] = i == j ? A[i][i] : 0.0;
			Drev[i][j] = i == j ? 1 / A[i][i] : 0.0;
			DminusA[i][j] = D[i][j] - A[i][j];
		}
	}


	double** B = productOfMatricies(N, Drev, DminusA);

	double* g = productOfMatrixAndVector(N, Drev, b);

	double* xk = new double[N];
	for (int i = 0; i < N; i++)
	{
		xk[i] = b[i];
	}
	double* xk1 = new double[N];
	const double eps = 1e-18;
	double q = 0.0;
	do
	{
		double* Bxk = productOfMatrixAndVector(N, B, xk);

		for (int i = 0; i < N; i++)
		{
			xk1[i] = Bxk[i] + g[i];
		}

		xk = xk1;

		double* Axk = productOfMatrixAndVector(N, A, xk);
		q = 0.0;
		for (int i = 0; i < N; i++)
		{
			q += (Axk[i] - b[i]) * (Axk[i] - b[i]);
		}

	} while (sqrt(q) > eps);

	for (int i = 0; i < N; i++)
	{
		x[i] = xk[i];
	}
}

/**
 * Метод минимальных невязок
 */
void nuyanzinma::lab6()
{
	double* xk = b;
	double* rk = new double[N];
	double tau = 0;
	double eps = 1e-18;
	double maxDiff;
	do
	{
		double* product = productOfMatrixAndVector(N, A, xk);
		for (int i = 0; i < N; i++)
		{
			rk[i] = product[i] - b[i];
		}

		//Скалярное произведение
		double* Ark = productOfMatrixAndVector(N, A, rk);
		double scalarProd1 = 0;
		for (int i = 0; i < N; i++)
		{
			scalarProd1 += Ark[i] * rk[i];
		}

		double scalarProd2 = 0;
		for (int i = 0; i < N; i++)
		{
			scalarProd2 += Ark[i] * Ark[i];
		}

		tau = scalarProd1 / scalarProd2;

		// Поиск Xk+1
		double* xk1 = new double[N];
		for (int i = 0; i < N; i++)
		{
			xk1[i] = xk[i] - tau * rk[i];
		}
		maxDiff = abs(xk1[0] - xk[0]);

		for (int i = 1; i < N; i++)
		{
			double currDiff = abs(xk1[i] - xk[i]);
			if (currDiff > maxDiff)
			{
				maxDiff = currDiff;
			}
		}
		xk = xk1;
	} while (maxDiff > eps);

	for (int i = 0; i < N; i++)
	{
		x[i] = xk[i];
	}
}

double scalarProduct(int N, double*& a, double*& b)
{
	double result = 0;
	for (int i = 0; i < N; i++)
	{
		result += a[i] * b[i];
	}

	return result;
}

double* sumOfVectors(int N, double*& a, double*& b, double a1, double b1)
{
	double* result = new double[N];

	for (int i = 0; i < N; i++)
	{
		result[i] = a1 * a[i] + b1 * b[i];
	}

	return result;
}

/**
 * Метод сопряженных градиентов
 */
void nuyanzinma::lab7()
{
	const double eps = 1e-40;
	double* buffer = productOfMatrixAndVector(N, A, b); // Вычисление Ax0, изначально x0 = b
	double* rk = new double[N];
	for (int i = 0; i < N; i++) // Вычисление r0
	{
		rk[i] = b[i] - buffer[i];
	}
	double* zk = rk; // Вычисление z0
	double* xk = b;
	double* xk1 = new double[N];
	double* rk1 = new double[N];
	double* zk1 = new double[N];
	double otnNev;

	do 
	{
		buffer = productOfMatrixAndVector(N, A, zk);
		double alphaK = scalarProduct(N, rk, rk) / scalarProduct(N, buffer, zk);
		xk1 = sumOfVectors(N, xk, zk, 1, alphaK);
		rk1 = sumOfVectors(N, rk, buffer, 1, -alphaK);
		double betaK = scalarProduct(N, rk1, rk1) / scalarProduct(N, rk, rk);
		zk1 = sumOfVectors(N, rk1, zk, 1, betaK);
		
		zk = zk1;
		rk = rk1;
		xk = xk1;
		
		otnNev = sqrt(scalarProduct(N, rk, rk) / scalarProduct(N, b, b));
	} while (abs(otnNev) > eps);

	for (int i = 0; i < N; i++)
	{
		x[i] = xk[i];
	}
}

double** transposeMatrix(double N, double**& M)
{
	double** tM = new double*[N];
	for (int i = 0; i < N; i++)
	{
		tM[i] = new double[N];
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			tM[i][j] = M[j][i];
		}
	}

	return tM;
}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void nuyanzinma::lab8()
{
	double eps = 1e-1;

	double** U = new double* [N];
	for (int i = 0; i < N; i++)
	{
		U[i] = new double[N];
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			U[i][j] = i == j ? 1 : 0;
		}
	}

	bool flag = false;
	do
	{
		flag = false;
		int iF = 0, jF = 0;

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				if (i != j && abs(A[i][j]) > abs(A[iF][jF]))
				{
					iF = i; jF = j;
				}
			}
		}

		double p = 2 * A[iF][jF] / (A[iF][iF] - A[jF][jF]);
		double s = sin(0.5 * atan(p));
		double c = cos(0.5 * atan(p));

		U[iF][iF] = c; U[iF][jF] = -s;
		U[jF][iF] = s; U[jF][jF] = c;

		double** trU = transposeMatrix(N, U);

		double** trUonA = productOfMatricies(N, trU, A);

		A = productOfMatricies(N, trUonA, U);
	
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				if (i != j && abs(A[i][j]) > eps)
				{
					flag = true;
				}
			}
		}

		U[iF][iF] = 1; U[iF][jF] = 0;
		U[jF][iF] = 0; U[jF][jF] = 1;

	} while (flag);

	for (int i = 0; i < N; i++)
	{
		x[i] = A[i][i];
	}
}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void nuyanzinma::lab9()
{
	double EPS = 1e-3;
	double* prevY = new double[N];
	for (int i = 0; i < N; i++)
	{
		prevY[i] = 1;
	}
	double* nextY = new double[N];
	double D = 0;
	double prevLambda = 0;
	double nextLambda = 0;

	do {
		nextY = productOfMatrixAndVector(N, A, prevY);

		nextLambda = nextY[N - 1] / prevY[N - 1];

		D = abs(nextLambda - prevLambda);

		prevLambda = nextLambda;
		prevY = nextY;

	} while (D > EPS);

	cout << nextLambda << endl;
}


std::string nuyanzinma::get_name()
{
  return "M.A. Nuyanzin";
}
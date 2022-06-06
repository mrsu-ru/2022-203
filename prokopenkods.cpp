#include "prokopenkods.h"

/**
 * Введение в дисциплину
 */
void prokopenkods::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void prokopenkods::lab2()
{
	for (int i = 0; i < N; i++) {

		int maxNum = i;

		for (int j = i + 1; j < N; j++) {
			if (abs(A[j][i]) > abs(A[maxNum][i])) {
				maxNum = j;
			}
		}

		if (maxNum != i) {
			for (int j = 0; j < N; j++) {
				swap(A[maxNum][j], A[i][j]);
				swap(b[maxNum], b[i]);
			}
		}

		for (int j = i + 1; j < N; j++) {
			double s = A[j][i] / A[i][i];
			for (int c = i; c < N; c++) {
				A[j][c] -= s * A[i][c];
			}
			b[j] -= (s * b[i]);
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
void prokopenkods::lab3()
{
	double* alpha = new double[N];
	double* beta = new double[N];
	double y;

	y = A[0][0];
	alpha[0] = -A[0][1] / y;
	beta[0] = b[0] / y;

	for (int i = 1; i < N; i++) {
		y = A[i][i] + A[i][i - 1] * alpha[i - 1];
		alpha[i] = -A[i][i + 1] / y;
		beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / y;
	}
	x[N - 1] = beta[N - 1];

	for (int i = N - 2; i >= 0; i--) {
		x[i] = alpha[i] * x[i + 1] + beta[i];
	}
}



/**
 * Метод Холецкого
 */
void prokopenkods::lab4()
{
	double S[N][N];
	int D[N][N];
	for (int i = 0; i < N; i++) {

		double subSum = 0;
		for (int l = 0; l <= i - 1; l++) {
			subSum += S[l][i] * S[l][i] * D[l][l];
		}

		int sign = (A[i][i] - subSum) < 0;
		D[i][i] = (int)pow(-1, sign);

		S[i][i] = sqrt(abs(A[i][i] - subSum));
		for (int j = 0; j < i; j++) S[i][j] = 0.0;
		for (int j = i + 1; j < N; j++) {
			double subSum = 0;
			for (int l = 0; l <= i - 1; l++) subSum += S[l][j] * S[l][i] * D[l][l];
			S[i][j] = (A[i][j] - subSum) / S[i][i] * D[i][i];
		}

		//Обратный ход

		double y[N];
		y[0] = b[0] / S[0][0];
		for (int i = 1; i < N; i++) {
			double subSum = 0;
			for (int j = 0; j <= i - 1; j++) subSum += S[j][i] * y[j];
			y[i] = b[i] - subSum;
			y[i] /= S[i][i];
		}
		x[N - 1] = y[N - 1] / S[N - 1][N - 1];
		for (int i = N - 2; i >= 0; i--) {
			double subSum = 0;
			for (int k = i + 1; k <= N - 1; k++) subSum += S[i][k] * x[k];
			x[i] = y[i] - subSum;
			x[i] /= S[i][i];
		}

	}

}



/**
 * Метод Якоби или Зейделя
 */
void prokopenkods::lab5()
{
	double* xk = new double[N];
	double eps = 1e-20;
	double norm;
	for (int i = 0; i < N; i++) x[i] = b[i] / A[i][i];

	do {
		for (int i = 0; i < N; i++) {

			for (int j = 0; j < N; j++) xk[j] = x[j];

			double lowerSum = 0, upperSum = 0;

			for (int j = 0; j < i; j++) lowerSum += A[i][j] * xk[j];
			for (int j = i + 1; j < N; j++) upperSum += A[i][j] * xk[j];

			x[i] = 1 / A[i][i] * (b[i] - lowerSum - upperSum);

			if (i == 0) norm = abs(x[i] - xk[i]);
			if (abs(x[i] - xk[i]) > norm) norm = abs(x[i] - xk[i]);

		}

	} while (norm >= eps);
}



/**
 * Метод минимальных невязок
 */
void prokopenkods::lab6()
{
	const double EPS = 1e-19;
	double z = 1e9;
	const int n = N;
	double r[n];
	int iter;
	for (iter = 0; z > EPS; iter++) {
		z = 0;
		for (int i = 0; i < n; i++) {
			r[i] = -b[i];
			for (int j = 0; j < n; j++) {
				r[i] += A[i][j] * x[j];
			}
		}
		double tlower = 0, tupper = 0;
		for (int i = 0; i < n; i++) {
			double temp = 0;
			for (int j = 0; j < n; j++) {
				temp += A[i][j] * r[j];
			}
			tlower += temp * temp;
			tupper += temp * r[i];
		}
		double t = tupper / tlower;
		for (int i = 0; i < n; i++) {
			if (z < abs(t * r[i])) z = abs(t * r[i]);
			x[i] -= t * r[i];
		}
	}
}



/**
 * Метод сопряженных градиентов
 */
void prokopenkods::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void prokopenkods::lab8()
{
	const double eps = 1e-15;
	double** C = new double* [N];
	for (int i = 0; i < N; i++) C[i] = new double[N];

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			C[i][j] = A[i][j];
		}
	}

	double norm = 1e9;
	int iter;
	for (iter = 0; norm > eps; iter++) {
		int k = 1, l = 2;
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				if (abs(A[k][l]) < abs(A[i][j])) {
					k = i;
					l = j;
				}
			}
		}

		double phi;
		if (fabs(A[k][k]-A[l][l]) < eps) {
			phi = atan(1);
		}
		else {
			phi = 0.5*atan(2*A[k][l]/(A[l][l]-A[k][k]));
		}

		double s = sin(phi);
		double c = cos(phi);

		C[k][k] = c*c*A[k][k]-2*s*c*A[k][l]+s*s*A[l][l]; // cos^2-2sincos+sin^2 = ajj
		C[l][l] = s*s*A[k][k]+2*s*c*A[k][l]+c*c*A[l][l]; // sin^2+2sincos+cos^2 = aii
		C[k][l] = C[l][k] = (c*c-s*s)*A[k][l]+s*c*(A[k][k]-A[l][l]); // cos^2-sin^2+sincos = aji


		for (int i = 0; i < N; i++) {
			if (i == k || i == l) continue;

			C[k][i] = C[i][k] = c*A[k][i]-s*A[l][i];
			C[l][i] = C[i][l] = s*A[k][i]+c*A[l][i];
		}

		norm = 0;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				A[i][j] = C[i][j];
				if (i < j) norm += A[i][j] * A[i][j];
			}
		}
	}


	for (int i = 0; i < N; i++) {
		x[i] = A[i][i];
	}

	for (int i = 0; i < N; i++) {
		delete[] C[i];
	}
	delete[] C;
}

/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void prokopenkods::lab9()
{
	const double eps = 1e-3;

	double* yPrev = new double[N];
	for (int i = 0; i < N; i++)
		yPrev[i] = 1;

	double* y = new double[N];

	int iter;
	double delta = 1e9;
	double lambdaMax = 0;

	for (iter = 0; delta > eps; iter++) {
		for (int i = 0; i < N; i++) {
			y[i] = 0;
			for (int j = 0; j < N; j++) {
				y[i] += A[i][j] * yPrev[j];
			}
		}

		double sum1 = 0, sum2 = 0;
		for (int i = 0; i < N; ++i) {
			sum1 += y[i];
			sum2 += yPrev[i];
		}
		double lambda = sum1 / sum2;
		delta = fabs(lambda - lambdaMax);
		lambdaMax = lambda;

		for (int i = 0; i < N; i++) {
			yPrev[i] = y[i];
		}
		double r = 0;
		for (int i = 0; i < N; i++) {
			r += y[i] * y[i];
		}
		r = sqrt(r);
		for (int i = 0; i < N; ++i) {
			y[i] /= r;
		}
	}
	printf("Max Eigen value -- %.4f\n", lambdaMax);

	delete[]yPrev;
	delete[]y;
}


std::string prokopenkods::get_name()
{
  return "D.S.Prokopenko";
}

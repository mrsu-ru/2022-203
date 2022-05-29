#include "kochetkovpa.h"
#include <iostream>
#include <math.h>

#define eps 1e-20

/**
 * Введение в дисциплину
 */
void kochetkovpa::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void kochetkovpa::lab2()
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
void kochetkovpa::lab3()
{
	double* alpha = new double[N];
	double* beta = new double[N];
	double y;
	
	y = A[0][0];
	alpha[0] = -A[0][1] / y;
	beta[0] = b[0] / y;
	
	for (int i = 1; i < N; i++) {
		y = A[i][i] + A[i][i - 1] * alpha[i - 1];
		alpha[i] = -A[i][i+1] / y;
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
void kochetkovpa::lab4()
{
    double S[N][N];
    int D[N][N];
    for (int i = 0; i < N; i++) {

        double subSum = 0;
        for (int l = 0; l <= i - 1; l++) {
            subSum += S[l][i] * S[l][i] * D[l][l];
        }

        int sign = (A[i][i] - subSum) < 0;
        D[i][i] =  (int)pow(-1, sign);

        S[i][i] = sqrt(abs(A[i][i] - subSum));
        for (int j = 0; j < i; j++) S[i][j] = 0.0;
        for (int j = i+1; j < N; j++) {
            double subSum = 0;
            for (int l = 0; l <= i - 1; l++) subSum += S[l][j] * S[l][i] * D[l][l];
            S[i][j] = (A[i][j] - subSum)/S[i][i]*D[i][i];
        }

        double y[N];
        y[0] = b[0] / S[0][0];
        for (int i = 1; i < N; i++) {
              double subSum = 0;
              for (int j = 0; j <= i - 1; j++) subSum += S[j][i] * y[j];
              y[i] = b[i] - subSum;
              y[i] /= S[i][i];
        }
        x[N - 1] = y[N - 1] / S[N - 1][N - 1];
        for (int i = N - 2; i >=0; i--) {
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
void kochetkovpa::lab5()
{
    double xk[N];
    double norm;
    double LU[N][N];

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (i != j)
            {
                LU[i][j] = A[i][j];
            }

            if (i == j)
            {
                LU[i][j] = 0.;
            }
        }
    }

    for (int i = 0; i < N; ++i)
        x[i] = b[i] / A[i][i];

    do
    {
        for (int i = 0; i < N; ++i)
        {
            for (int i = 0; i < N; ++i)
                xk[i] = x[i];

            double summ = 0.;

            for (int j = 0; j < N; ++j)
            {
                summ += LU[i][j] * xk[j];
            }

            x[i] = (1 / A[i][i]) * (b[i] - summ);

            if (i == 0)
            {
                norm = fabs(x[i] - xk[i]);
            }

            if (fabs(x[i] - xk[i]) > norm)
            {
                norm = fabs(x[i] - xk[i]);
            }
        }
    } while (sqrt(norm) >= eps);
}


/**
 * Метод минимальных невязок
 */
void kochetkovpa::lab6()
{
    double norm = 1e-1;
    double rk[N];
    init_matrix(rk, N);

    for (int i = 0; i < N; ++i)
    {
        x[i] = 1e-2;
    }

    do
    {
        double *Ax = MatrixMultOnVect(A, x, N);

        for (int i = 0; i < N; ++i)
        {
            rk[i] = Ax[i] - b[i];
        }

        double *Ar = MatrixMultOnVect(A, rk, N);
        double tau = 0;

        tau = (scalar(Ar, rk, N) / scalar(Ar, Ar, N));

        for (int i = 0; i < N; i++)
        {
            double prevX = x[i];
            x[i] = prevX - tau * rk[i];

            if (fabs(x[i] - prevX) < norm)
            {
                norm = fabs(x[i] - prevX);
            }
        }

        delete[] Ar;
        delete[] Ax;

    } while (sqrt(norm) >= eps);
}



/**
 * Метод сопряженных градиентов
 */
void kochetkovpa::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void kochetkovpa::lab8()
{
  double err = 0;
  for (int i = 0; i < N; i++) 
    for (int j = i + 1; j < N; j++)
      err += A[i][j] * A[i][j] + A[j][i] * A[j][i];
  
  double **C = new double*[N];
  for (int i = 0; i < N; i++)
    C[i] = new double[N];

  while (err > eps)
  {
    double alpha = 0;

    int i = -1, j = -1;
    double max = -1e9;
    for (int ii = 0; ii < N; ii++)
    {
      for (int jj = ii + 1; jj < N; jj++)
      {
        if (fabs(A[ii][jj]) > max)
        {
           max = fabs(A[ii][jj]);
           i = ii; j = jj; 
        }
      }
    }

    if (fabs(A[i][i] - A[j][j]) < eps)   
      alpha = atan(1);
    else
      alpha = atan(2 * A[i][j] / (A[j][j] - A[i][i])) / 2;

    double s = sin(alpha), c = cos(alpha);

    for (int ii = 0; ii < N; ii++)
      for (int jj = 0; jj < N; jj++)
          C[ii][jj] = A[ii][jj];     
      
    C[i][i] = c * c * A[i][i] - 2 * s * c * A[i][j] + s * s * A[j][j];
    C[j][j] = s * s * A[i][i] + 2 * s * c * A[i][j] + c * c * A[j][j];
    C[i][j] = C[j][i] = (c * c - s * s) * A[i][j] + s * c * (A[i][i] - A[j][j]);
    for (int k = 0; k < N; k++)
    {
      if (k == i || k == j) continue;
      C[i][k] = C[k][i] = c * A[i][k] - s * A[j][k];
      C[j][k] = C[k][j] = s * A[i][k] + c * A[j][k];
    }

    err = 0;
    for (int i = 0; i < N; i++) 
      for (int j = i + 1; j < N; j++)
        err += C[i][j] * C[i][j] + C[j][i] * C[j][i];

    for (int ii = 0; ii < N; ii++)
      for (int jj = 0; jj < N; jj++)
          A[ii][jj] = C[ii][jj];
  }

  for (int i = 0; i < N; i++)
  {
      cout  << "lambda " << i << " = " << A[i][i] << endl;  
      delete []C[i];
  }
  
  delete []C;
}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void kochetkovpa::lab9()
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            if (A[i][j] != A[j][i]) return;
    }
    double l, maxl = 0;
    bool flag = true;
    double *y = new double[N];
    double *y_prev = new double[N];
    for (int i = 0; i < N; i++) y_prev[i] = 1;
    while (flag) {
        flag = false;
        for (int i = 0; i < N; i++) {
            y[i] = 0;
            for (int j = 0; j < N; j++) {
                y[i] += A[i][j] * y_prev[j];
            }
        }
        l = 0;
        for (int i = 0; i < N; i++) {
            if (fabs(y[i]) > eps && fabs(y_prev[i]) > eps) {
                l = y[i] / y_prev[i];
                break;
            }
        }
        if (fabs(l - maxl) > eps) flag = true;
        maxl = l;
        for (int i = 0; i < N; i++) y_prev[i] = y[i];
    }
    cout << "The largest lambda: " << maxl << endl;
    delete[] y;
    delete[] y_prev;
}


std::string kochetkovpa::get_name()
{
  return "P.A. Kochetkov";
}

#include "fedinda.h"
#include <iostream>
const double eps = 1.e-15;
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
	double* alpha = new double[N];
	double* beta = new double[N];
	int i;
	alpha[0] = A[0][1] / A[0][0];
	beta[0] = b[0] / A[0][0];
	for (i = 1; i < N; i++) {
		alpha[i] = A[i][i + 1] / (A[i][i] - A[i][i - 1] * alpha[i - 1]);
		beta[i] = (b[i] - beta[i - 1] * A[1][0]) / (A[i][i] - A[i][i - 1] * alpha[i - 1]);
	}

	x[N - 1] = beta[N - 1];
	for (i = N - 2; i >= 0; i--) {
		x[i] = beta[i] - alpha[i] * x[i + 1];
	}
}



/**
 * Метод Холецкого
 */
void fedinda::lab4()
{
	const double lambda = 0.01;
	double prev[N];
	for (int i = 0; i < N; i++){
		prev[i] = 0;
	}

	double difference;
	do{
		for (int i = 0; i < N; i++){
			double sum = 0;
			for (int j = 0; j < N; j++){
				sum += A[i][j] * prev[j];
			}
			x[i] = prev[i] - lambda * (sum - b[i]);
		}
		difference = 0;
		for (int i = 0; i < N; i++){
			if (fabs(x[i] - prev[i]) > difference){
				difference = fabs(x[i] - prev[i]);
			}
		}
		for (int i = 0; i < N; i++){
			prev[i] = x[i];
		}
	} while (difference > eps);
}



/**
 * Метод Якоби или Зейделя
 */
void fedinda::lab5()
{
	double prev[N];
	for (int i = 0; i < N; i++)
	{
		prev[i] = 0;
	}

	double difference;
	do
	{
		for (int i = 0; i < N; i++)
		{
			double sum = 0;
			for (int j = 0; j < i; j++)
			{
				sum += A[i][j] * prev[j];
			}

			for (int j = i + 1; j < N; j++)
			{
				sum += A[i][j] * prev[j];
			}

			x[i] = (b[i] - sum) / A[i][i];
		}

		difference = 0;
		for (int i = 0; i < N; i++)
		{
			if (fabs(x[i] - prev[i]) > difference)
			{
				difference = fabs(x[i] - prev[i]);
			}
		}

		for (int i = 0; i < N; i++)
		{
			prev[i] = x[i];
		}

	} while (difference > eps);
}



/**
 * Метод минимальных невязок
 */
void fedinda::lab6()
{
    double err;
    double Ax, Ay[N];
    double y0[N];
    double t, x_prev;
    double sum1, sum2;
    do{
        for (int i=0; i<N; i++){
            Ax=0;
            for (int j=0; j<N; j++){
                Ax+=A[i][j]*x[j];
            }
            y0[i]=b[i]-Ax;
        }
        for (int i=0; i<N; i++){
            Ay[i]=0;
            for (int j=0; j<N; j++){
                Ay[i]+=A[i][j]*y0[j];
            }
        }
        sum1=0; sum2=0;
        for (int i=0; i<N; i++){
            sum1+=y0[i]*Ay[i];
            sum2+=Ay[i]*Ay[i];
        }
        t=sum1/sum2;
        err=0;
        for(int i=0; i<N; i++){
            x_prev=x[i];
            x[i]+=t*y0[i];
            if (abs(x[i]-x_prev)>err){
                err=abs(x[i]-x_prev);
            }
        }
    }while (err>eps);
}



/**
 * Метод сопряженных градиентов
 */
void fedinda::lab7()
{
	double *r = new double[N];
	double *z = new double[N];
	double *Az = new double[N];
	double alfa, betta, tmp, norma, b_norm=0;

	for(int i=0; i<N; i++){
		b_norm += b[i]*b[i]; 
		r[i] = b[i];
		for(int j=0; j<N; j++)
			r[i]-=A[i][j]*x[j];
		z[i]=r[i];
	}

	do{
		for(int i=0; i<N; i++){
			Az[i]=0;
			for(int j=0; j<N; j++)
				Az[i]+=A[i][j]*z[j];
		}
		alfa=0; tmp=0;
		for(int i=0; i<N; i++){
			alfa+=r[i]*r[i];
			tmp+=Az[i]*r[i];
		}
		alfa/=tmp;

		betta=0; tmp=0;
    	for(int i=0; i<N; i++){
			tmp+=r[i]*r[i];
			x[i]+=alfa*z[i];
			r[i]-=alfa*Az[i];
			betta+=r[i]*r[i];
		}
		norma=sqrt(betta/b_norm);
		betta/=tmp;
		
		for(int i=0; i<N; i++)
			z[i]=r[i]+betta*z[i];
	}while(norma>eps);
}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void fedinda::lab8()
{
	double t = 1;
	double** B = new double* [N];
	B[0] = new double[N * N];
	for (int i = 0; i < N; i++)
		B[i] = B[0] + i * N;

	while (t > eps) {
		int i_fix = 1;
		int j_fix = 0;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < i; j++)
				if (abs(A[i][j]) > abs(A[i_fix][j_fix])) {
					i_fix = i;
					j_fix = j;
				}

		double phi = 0.5 * atan(2 * A[i_fix][j_fix] / (A[i_fix][i_fix] - A[j_fix][j_fix]));

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				B[i][j] = A[i][j];

		for (int i = 0; i < N; i++)
		{
			B[i][i_fix] = A[i][i_fix] * cos(phi) + A[i][j_fix] * sin(phi);
			B[i][j_fix] = A[i][j_fix] * cos(phi) - A[i][i_fix] * sin(phi);
		}
		for (int i = 0; i < N; i++)
		{
			A[i][i_fix] = B[i][i_fix];
			A[i][j_fix] = B[i][j_fix];
		}
		for (int j = 0; j < N; j++)
		{
			A[i_fix][j] = B[i_fix][j] * cos(phi) + B[j_fix][j] * sin(phi);
			A[j_fix][j] = B[j_fix][j] * cos(phi) - B[i_fix][j] * sin(phi);
		}
		t = 0;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < i; j++)
				t += A[i][j] * A[i][j];
		t *= 2;
	}
	for (int i = 0; i < N; i++)
		x[i] = A[i][i];

	delete[] B[0];
	delete[] B;
}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void fedinda::lab9()
{
	double* y = new double[N];
	double* prev_y = new double[N];
	double max_lambda = 0;
	for (int i = 0; i < N; i++)
		prev_y[i] = 1;
	bool stop = false;

	while (!stop)
	{
		for (int i = 0; i < N; i++)
		{
			y[i] = 0;
			for (int j = 0; j < N; j++)
				y[i] += A[i][j] * prev_y[j];
		}
		double tmp = max_lambda;
		for (int i = 0; i < N; i++)
		{
			if (abs(y[i]) > eps && abs(prev_y[i]) > eps)
			{
				max_lambda = y[i] / prev_y[i];
				break;
			}
		}

		for (int i = 0; i < N; i++)
			prev_y[i] = y[i];

		if (abs(max_lambda - tmp) < eps)
			stop = true;

	}
	cout << "max eigenvalue: " << max_lambda;

	delete[] prev_y;
	delete[] y;
}


std::string fedinda::get_name()
{
  return "D.A. Fedin";
}

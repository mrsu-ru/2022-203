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
	double err = 0;
	double** C = new double*[N];

	for (int i = 0; i < N; i++){
		C[i] = new double[N];
	}


	for (int i = 0; i < N; i++){
		for (int j = i+1; j < N; j++){
    			if(i!=j){ 
      				err+=A[i][j]* A[i][j];
   			}
    			C[i][j] = 0;
  		}
	}

	while (sqrt(err) > eps){
		int mI = 0, mJ = 1;
		for (int i = 0; i < N; i++)
			for (int j = i+1; j<N; j++)
				if ( abs(A[i][j]) > abs(A[mI][mJ]) ){
  	    				mI = i; mJ = j;	
				}

	double phi;
  	if(A[mI][mI]!= A[mJ][mJ]){
    		phi = 0.5*atan(2*A[mI][mJ]/(A[mI][mI] - A[mJ][mJ]));
 	}else{phi = M_PI/4;}
	double c = cos(phi), s = sin(phi);
 
	C[mI][mI] = pow(c, 2)*A[mI][mI] - 2*s*c*A[mI][mJ] + pow(s, 2)*A[mJ][mJ];
	C[mJ][mJ] = pow(s, 2)*A[mI][mI] + 2*s*c*A[mI][mJ] + pow(c, 2)*A[mJ][mJ];
	C[mI][mJ] = (pow(c, 2) - pow(s, 2))*A[mI][mJ] + s*c*(A[mJ][mJ] - A[mI][mI]);
	C[mJ][mI] = C[mI][mJ];
	
	for (int k = 0; k < N; k++){
		if (k != mI  &&  k != mJ){
	  		C[mI][k] = c*A[mI][k] - s*c*A[mJ][k];
	  		C[k][mI] = C[mI][k];
	  		C[mJ][k] = s*A[mI][k] + c*A[mJ][k];
	  		C[k][mJ] = C[mJ][k];
		} 
		for (int l = 0; l < N; l++)
	    		if (k != mI && k != mJ && l != mI && l != mJ) 
				C[k][l] = A[k][l];	  
	}
	
	err = 0;
	for (int i = 0; i < N; i++)
        	for (int j = i+1; j < N; j++)
  	        	if (i != j) 
			err += C[i][j] * C[i][j]; 

	for (int i = 0; i < N; i++)
        	for (int j = 0; j < N; j++) 
	  		A[i][j] = C[i][j]; 

  }
  
  	for (int i = 0; i < N; i++) 
  		x[i] = A[i][i];
}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void fedinda::lab9()
{
	int n = N;
	double *y = new double[n];
	double *y_next = new double[n];
	double lyambda = 1;
  	double lyambda_next = 0;

  	for(int i = 0; i<n; i++)
    		y[i] = b[i];

  	do{
  		for(int i=0; i<n; i++){
      			for(int j=0; j<n; j++){
        			y_next[i] += A[i][j]*y[j];
      			}
  		}
  		lyambda = lyambda_next;

  		for(int i=0; i<n; i++){
    			if(y[i]!= 0 && y_next[i] != 0){
      				lyambda_next = y_next[i]/y[i];
      				break;
    			}
  		}
  		for (int i=0; i<n; i++)
   			y[i]=y_next[i];
	}while(fabs(lyambda_next - lyambda)>eps);

	cout<<"Result: "<<lyambda_next << endl;
}


std::string fedinda::get_name()
{
  return "D.A. Fedin";
}

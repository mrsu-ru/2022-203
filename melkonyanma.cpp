#include "melkonyanma.h"

void MatrixTrans(double** matrix,double** tMatrix, int N);
void init_matrix(double* matrix, int N);
void init_matrix(double** matrix, int N);
void create_dynamic_array(double** array,int N);
void delete_dynamic_array(double** array,int N);

/**
 * Введение в дисциплину
 */
void melkonyanma::lab1()
{
  cout << "Hello World!" << endl;
}

void method_Gauss(double** A,double* b,int N)
{
    for(int row = 0, column = 0 ; column < N ; ++row, ++column)
    {
        int max = row;

        for(int i = row + 1; i < N ;++i)
        {
            if(fabs(A[i][column]) > fabs(A[max][column]))
            {
                max = i;
            }
        }

        if(max != row){
            swap(A[row],A[max]);
            swap(b[row],b[max]);
        }

        double mainElem = A[row][column];

        for(int i = 0; i < N ;++i)
        {
            if(A[row][i] == 0){ continue; }
            A[row][i] /= mainElem;           
        }

        b[row] /= mainElem; 

        for(int i = row + 1; i < N ;++i)
        {
            double a = A[i][column];
            for(int j = 0; j < N ;++j){
                A[i][j] -= A[row][j]*a;
            }
            b[i] -=  b[row]*a;
        }
    }   

    for(int i = N - 1; i > 0 ; --i)
    {
        double a = A[i][i];
        for(int k = i-1; k >= 0 ;--k)
        {
            b[k] -= b[i]*A[k][i];
            A[k][i] -= a*A[k][i];
        }
    } 
}

/**
 * Метод Гаусса с выбором главного элемента
 */
void melkonyanma::lab2()
{
    method_Gauss(A,b,N);

    for(int i = 0; i < N ; ++i)
    {
        x[i] = b[i];
    }
}



/**
 * Метод прогонки
 */
void melkonyanma::lab3()
{
    double * alpha = new double [N];
    double * beta = new double [N];
    double kappa;

    init_matrix(alpha,N);
    init_matrix(beta,N);

    kappa = A[0][0];
    alpha[0] = -A[0][1]/kappa;
    beta[0] = b[0]/kappa;

    for(int i = 1; i < N ; ++i)
    {
        kappa = A[i][i] + A[i][i - 1]*alpha[i - 1];
        
        if( i != N - 1 )
        {
            alpha[i] = -A[i][i + 1]/kappa;
        }

        beta[i] = (b[i] - A[i][i - 1]*beta[i - 1])/kappa;
    }

    x[N-1] = beta[N-1];

    for(int i = N - 2; i >= 0 ; --i)
    {
        x[i] = alpha[i]*x[i + 1] + beta[i];
    }
    
    delete [] alpha;
    delete [] beta;
}

/**
 * Метод Холецкого
 */
void melkonyanma::lab4()
{
    double* y = new double[N];
    double** L = new double*[N];
    double** Lt = new double*[N];
    double summ;

    create_dynamic_array(L,N);
    create_dynamic_array(Lt,N);

    init_matrix(L,N);
    init_matrix(Lt,N);
    init_matrix(y,N);

    L[0][0] = sqrt(A[0][0]);
    for(int i = 1; i < N;++i)
    {
        L[i][0] = A[i][0]/L[0][0];
    }

    for(int i = 1; i < N;++i)
    {   
        summ = 0.;

        for(int p = 0; p < i ;++p)
        {
            summ += L[i][p]*L[i][p]; 
        }
        L[i][i] = sqrt(A[i][i] - summ);

        summ = 0.;
        
        for(int j = i + 1; j < N ;++j)
        {
            for(int p = 0; p < i; ++p)
            {
                summ += L[i][p]*L[j][p];
            }

            L[j][i] = (A[j][i]-summ)/L[i][i];
        }    
    }

    MatrixTrans(L,Lt,N);
    method_Gauss(L,b,N);
    for(int i = 0; i < N ; ++i)
    {
        y[i] = b[i];
    }
    method_Gauss(Lt,y,N);
    for(int i = 0; i < N ; ++i)
    {
        x[i] = y[i];
    }

   
    delete_dynamic_array(L,N);
    delete_dynamic_array(Lt,N);

    delete [] y;
}


/**
 * Метод Якоби или Зейделя
 */
void melkonyanma::lab5()
{

}



/**
 * Метод минимальных невязок
 */
void melkonyanma::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void melkonyanma::lab7()
{

}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void melkonyanma::lab8()
{

}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void melkonyanma::lab9()
{

}

void MatrixTrans(double** matrix,double** tMatrix, int N)
{
    for(int i = 0; i < N ; ++i)
    {
        for(int j = 0; j < N ;++j)
        {
            tMatrix[i][j] = matrix[j][i];
        }
    }
}

void init_matrix(double* matrix, int N)
{
     for(int i = 0; i < N ;++i)
    {
            matrix[i] = 0.;
    }
}

void init_matrix(double** matrix,int N)
{
    for(int i = 0; i < N ;++i)
    {
        for(int j = 0; j < N ;++j)
        {
            matrix[i][j] = 0.;
        }
    }
}

void create_dynamic_array(double** array,int N)
{
    for(int i = 0 ; i < N; ++i)
    {
        array[i] = new double[N];
    }
}

void delete_dynamic_array(double** array, int N)
{
    for(int i = 0 ; i < N; ++i)
    {
        delete [] array[i];
    }
}

std::string melkonyanma::get_name()
{
  return "M.A. Melkonyan";
}
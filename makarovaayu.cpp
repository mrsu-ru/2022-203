#include "makarovaayu.h"
#include "melkonyanma.h"

/**
 * Введение в дисциплину
 */
void makarovaayu::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void makarovaayu::lab2()
{
    for (int i = 0; i < N; i++) {
        int maxLine = i;
        for (int j = i + 1; j < N; j++) {
            if (abs(A[j][i]) > abs(A[maxLine][i])) {
                maxLine = j;
            }
        }
        if (i != maxLine) {
            swap(A[maxLine], A[i]);
            swap(b[maxLine], b[i]);
        }
        for (int j = i + 1; j < N; j++) {
            double c = A[j][i] / A[i][i];
            for (int k = i; k < N; k++) {
                A[j][k] -= (c * A[i][k]);
            }
            b[j] -= (c * b[i]);
        }
    }

    for (int i = N - 1; i >= 0; i--) {
        for (int j = i + 1; j < N; j++) {
            b[i] -= (A[i][j] * x[j]);
        }
        x[i] = b[i] / A[i][i];
    }

}



/**
 * Метод прогонки
 */
void makarovaayu::lab3()
{
        double alpha[N - 1], beta[N];
        alpha[0] = - A[0][1] / A[0][0];
        beta[0] = b[0] / A[0][0];

        for (int i = 1; i < N; i++)
        {
            double y = A[i][i] + A[i][i - 1] * alpha[i - 1];
            alpha[i] = - A[i][i + 1] / y;
            beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / y;
        }

        x[N - 1] = beta[N - 1];
        for (int i = N - 2; i >= 0; i--)
        {
            x[i] = alpha[i] * x[i + 1] + beta[i];
        }
    }





/**
 * Метод Холецкого
 */
void makarovaayu::lab4()
{

        double **S = new double*[N];
        for (int i=0; i<N; i++)
        {
            S[i]=new double[N];
            for(int j=0; j<N; j++)
                S[i][j]=0;
        }
        double *D = new double[N];
        if (A[0][0]>0)
            D[0]=1;
        else
            D[0]=-1;
        S[0][0]=sqrt(fabs(A[0][0]));

        for (int i=1; i<N; i++)
        {
            S[0][i]=A[0][i]/(D[0]*S[0][0]);
        }

        for (int i=1; i<N; i++)
        {
            double temp =0;
            for (int j=0; j<i; j++)
                temp+=D[j]*S[j][i]*S[j][i];
            if (A[i][i]-temp>=0)
                D[i]=1;
            else
                D[i]=-1;
            S[i][i]=sqrt(D[i]*(A[i][i]-temp));

            for (int j=i+1; j<N; j++)
            {
                double l = 0;
                for (int k=0; k<j; k++)
                    l+=D[k]*S[k][i]*S[k][j];

                S[i][j]=(A[i][j]-l)/(D[i]*S[i][i]);
            }
        }
        double *y = new double[N];
        y[0]=b[0]/S[0][0];
        for (int i=1; i<N; i++)
        {
            double temp = 0;
            for (int j=0; j<i; j++)
                temp+=y[j]*S[j][i];
            y[i]=(b[i]-temp)/S[i][i];
        }
        x[N-1]=y[N-1]/(D[N-1]*S[N-1][N-1]);

        for (int i=N-2; i>=0; i--)
        {
            double temp =0;
            for (int j=i+1; j<N; j++)
                temp+=x[j]*D[j]*S[i][j];
            x[i]=(y[i]-temp)/(D[i]*S[i][i]);
        }



}



/**
 * Метод Якоби или Зейделя
 */
void makarovaayu::lab5()
{
    double norm, sum;
    double eps = 1.e-15;
    for (int i = 0; i < N; ++i) x[i] = b[i];

    do {
        norm = 0;
        for (int i = 0; i < N; i++) {
            sum = 0;
            for (int j = 0; j < N; j++) {
                if (i != j) sum += A[i][j] * x[j];
            }
            sum = (b[i] - sum) / A[i][i];
            if (norm < (fabs(sum - x[i]))) {
                norm = (fabs(sum - x[i]));
            }
            x[i] = sum;
        }
    } while (norm >= eps);
}


double *MulVecToMat(int N, double *A[], double b[]) {
    double *temp = new double[N];
    for (int i = 0; i < N; i++) {
        temp[i] = 0;
        for (int j = 0; j < N; j++) {
            temp[i] += A[i][j] * b[i];
        }
    }
    return temp;
}

double ScMul(int N, double temp[], double r[]) {
    double k = 0;
    for (int i = 0; i < N; i++) {
        k += temp[i] * r[i];
    }
    return k;
}



/**
 * Метод минимальных невязок
 */
void makarovaayu::lab6()
{
    double *r = new double[N];
    double *xk = new double[N];
    double eps = 1.e-17;
    double *x1 = b;
    double t;
    double maxDelta;
    do {
        double *temp = MulVecToMat(N, A, x1);

        for (int i = 0; i < N; i++) {
            r[i] = temp[i] - b[i];
        }
        double *Ar = MulVecToMat(N, A, r);

        double Scalar1, Scalar2;
        Scalar1 = ScMul(N, Ar, r);
        Scalar2 = ScMul(N, Ar, Ar);
        t = Scalar1 / Scalar2;

        for (int i = 0; i < N; i++) {
            xk[i] = x1[i] - t * r[i];
        }
        maxDelta = abs(xk[0] - x1[0]);

        for (int i = 1; i < N; i++) {
            double delta = abs(xk[i] - x1[i]);
            if (delta > maxDelta) {
                maxDelta = delta;
            }
        }
        x1 = xk;

    } while (maxDelta > eps);

    for (int i = 0; i < N; i++) {
        x[i] = x1[i];
    }
    delete[] r;
    delete[] xk;
}



/**
 * Метод сопряженных градиентов
 */
void makarovaayu::lab7() {
}
    double **TransposeMatrix(double **&m, int n) {
        double **temp = new double *[n];
        for (int i = 0; i < n; ++i) {
            temp[i] = new double[n];
        }

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                temp[i][j] = m[j][i];
            }
        }

        return temp;
    }

    double **MulMatrixes(double **&m1, double **&m2, int n) {
        double **temp = new double *[n];
        for (int i = 0; i < n; ++i) {
            temp[i] = new double[n];
        }

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                temp[i][j] = 0;
                for (int k = 0; k < n; ++k) {
                    temp[i][j] += m1[i][k] * m2[k][j];
                }
            }
        }

        return temp;

    }


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
    void makarovaayu::lab8() {
    double eps = 1.e-1;
    double MaxEl;

    do {
        MaxEl = 0;
        int maxi = 0, maxj = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i < j) {
                    if (MaxEl < fabs(A[i][j])) {
                        MaxEl = fabs(A[i][j]);
                        maxi = i;
                        maxj = j;
                    }
                }
            }
        }

        double Phi = 0.5 * atan(2 * A[maxi][maxj] / (A[maxi][maxi] - A[maxj][maxj]));

        double **m_H = new double *[N];
        for (int i = 0; i < N; ++i) {
            m_H[i] = new double[N];
        }

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i == j) {
                    m_H[i][j] = 1;
                } else {
                    m_H[i][j] = 0;
                }
            }
        }

        m_H[maxi][maxj] = -sin(Phi);
        m_H[maxj][maxi] = sin(Phi);

        m_H[maxi][maxi] = cos(Phi);
        m_H[maxj][maxj] = cos(Phi);

        double **t_H = TransposeMatrix(m_H, N);

        double **t_A = MulMatrixes(t_H, A, N);

        A = MulMatrixes(t_A, m_H, N);
    } while (MaxEl >= eps);

    for (int i = 0; i < N; ++i) {
        cout << "lambda  = " << A[i][i] << endl;
    }

    }


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
    void makarovaayu::lab9() {double AbsMaxEigenvalue;
    double Eigenvalue = 0.;
    double eps = 1e-8;
    double *xPrev = new double[N];
    double *xNew = new double[N];

    for (int i = 0; i < N; ++i)
    {
        xPrev[i] = b[i];
    }

    do
    {
        AbsMaxEigenvalue = Eigenvalue;

        xNew = MatrixMultVect(A, xPrev, N);

        double s1 = 0, s2 = 0;
        for (int i = 0; i < N; ++i)
        {
            s1 += xNew[i];
            s2 += xPrev[i];
        }

        Eigenvalue = s1 / s2;

        for (int i = 0; i < N; ++i)
        {
            xPrev[i] = xNew[i] / scalar(xNew, xNew, N);
        }

    } while (fabs(AbsMaxEigenvalue - Eigenvalue) >= eps);

    delete[] xNew;
    delete[] xPrev;

    cout << "Absolute max eigenvalue = " << AbsMaxEigenvalue << endl;
}




        std::string makarovaayu::get_name() {
            return "A.U.Makarova";
        }

double *makarovaayu::MatrixMultVect(double **pDouble, double *pDouble1, int n) {
    return nullptr;
}


#include "timovkinayu.h"

/**
 * Введение в дисциплину
 */
void timovkinayu::lab1()
{
  cout << "Hello World!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void timovkinayu::lab2()
{
    double p;
    int maxn;
    for (int k = 0; k < N - 1; k++) {
        maxn = k;
        for (int i = k + 1; i < N; i++)
            if (abs(A[i][k]) > abs(A[maxn][k])) maxn = i; 
        swap(A[maxn], A[k]);
        swap(b[maxn], b[k]);
        for (int i = k + 1; i < N; i++) {
            p = A[i][k] / A[k][k];
            for (int j = k; j < N; j++)
                A[i][j] -= p * A[k][j];
            b[i] -= p * b[k];
        }
    }
    for (int i = 0; i < N; i++) x[i] = b[i];
    for (int i = N - 1; i >= 0; i--) {
        for (int j = i + 1; j < N; j++)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
}



/**
 * Метод прогонки
 */
void timovkinayu::lab3()
{
	double* alpha = new double[N]; 
	double* betta = new double[N]; 
	alpha[0] = -A[0][1] / A[0][0];
	betta[0] = b[0] / A[0][0];
	for (int i = 1; i < N; i++) { 

		alpha[i] = A[i][i + 1] / (-A[i][i] - A[i][i - 1] * alpha[i - 1]);
		betta[i] = (-b[i] + A[i][i - 1] * betta[i - 1]) / (-A[i][i] - A[i][i - 1] * alpha[i - 1]);
	}
	x[N - 1] = betta[N - 1];
	for (int i = N - 2; i >= 0; i--) 
		x[i] = alpha[i] * x[i + 1] + betta[i];
	delete[] alpha;
	delete[] betta;
}



/**
 * Метод Холецкого
 */
void timovkinayu::lab4()
{
double *I = new double[N];
    double **W = new double*[N];
    for (int i = 0; i < N; i++) {
        W[i] = new double[N];
        for (int j = 0; j < N; j++) W[i][j] = 0;
    }
    if (A[0][0] > 0) I[0] = 1;
    else I[0] = -1;
    W[0][0] = sqrt(fabs(A[0][0]));
    for (int j = 1; j < N; j++) W[0][j] = A[0][j] / (I[0] * W[0][0]);
    for (int i = 1; i < N; i++) {
        double tmp = 0;
        for (int j = 0; j < i; j++) {
            tmp += I[j] * W[j][i] * W[j][i];
        }
        I[i] = copysign(1, A[i][i] - tmp);
        W[i][i] = sqrt(I[i] * (A[i][i] - tmp));

        for (int j = i + 1; j < N; j++) {
            double sum = 0;
            for (int k = 0; k < j; k++) {
                sum += I[k] * W[k][i] * W[k][j];
            }
            W[i][j] = (A[i][j] - sum) / (I[i] * W[i][i]);
        }
    }
    double* y = new double[N];
    y[0] = b[0] / W[0][0];
    for (int i = 1; i < N; i++) {
        double temp = 0;
        for (int j = 0; j < i; j++) temp += y[j] * W[j][i];
        y[i] = (b[i] - temp) / W[i][i];
    }
    x[N - 1] = y[N - 1] / (I[N - 1] * W[N - 1][N - 1]);
    for (int i = N - 2; i >= 0; i--) {
        double temp = 0;
        for (int j = i + 1; j < N; j++) temp += x[j] * I[j] * W[i][j];
        x[i] = (y[i] - temp) / (I[i] * W[i][i]);
    }
}



/**
 * Метод Якоби или Зейделя
 */
void timovkinayu::lab5()
{
double E = 1e-30;
    double* d = new double[N];
    double temp;
    do { 
        for (int i = 0; i < N; i++) d[i] = x[i];
        for (int i = 0; i < N; i++) {
            double sum1 = 0, sum2 = 0;
            for (int j = i + 1; j < N; j++) sum1 += A[i][j] * x[j];
            for (int j = i - 1; j >= 0; j--) sum2 += A[i][j] * x[j];
            x[i] = (b[i] - sum1 - sum2) / A[i][i]; }
        temp = 0;
        for (int i = 0; i < N; i++) temp += abs(x[i] - d[i]);
    } while (temp > E);
}



/**
 * Метод минимальных невязок
 */
void timovkinayu::lab6()
{
double* L = new double[N];
    double* R = new double[N];
    double nrm, eps = 1e-15;
    do {
        for (int i = 0; i < N; i++) {
            double temp = 0;
            for (int j = 0; j < N; j++) temp += A[i][j] * x[j];
            R[i] = temp - b[i];
            L[i] = 2 * R[i]; }
        double* A1 = new double[N];
        for (int i = 0; i < N; i++) {
            double temp = 0;
            for (int j = 0; j < N; j++) temp += A[i][j] * R[j];
            A1[i] = temp;
        }
        double t1 = 0, t2 = 0;
        for (int i = 0; i < N; i++) {
            t1 += abs(A1[i] * R[i]);
            t2 += abs(A1[i] * A1[i]); }
        double a = t1 / (2 * t2);
        double* y = new double[N];
        for (int i = 0; i < N; i++) y[i] = x[i];
        for (int i = 0; i < N; i++) x[i] = x[i] - a * L[i];
        nrm = 0;
        for (int i = 0; i < N; i++) nrm += (y[i] - x[i]) * (y[i] - x[i]);
    } while (sqrt(nrm) > eps);
}



/**
 * Метод сопряженных градиентов
 */
void timovkinayu::lab7()
{
int i, j; double e = 1.0e-21;
    double *Aq, *R, *q, *R1;
    q = new double[N]; Aq = new double[N];
    R = new double[N]; R1 = new double[N];
    for (int i = 0; i < N; i++) x[i] = b[i];
    for (i = 0; i < N; i++) {
        double s = 0;
        for (j = 0; j < N; j++) s += A[i][j] * x[j];
        R[i] = -s + b[i];
        q[i] = R[i]; }
    double nrm = 1; double k = 0;
    while (nrm > e) {
        double alfa1 = 0, alfa2 = 0;
        k++;
        for (i = 0; i < N; i++) alfa1 += R[i] * R[i];
        for (i = 0; i < N; i++) {
            double s = 0;
            for (j = 0; j < N; j++) s += A[i][j] * q[j];
            Aq[i] = s; }
        for (i = 0; i < N; i++) alfa2 += Aq[i] * q[i];
        double alfa = alfa1 / alfa2;
        for (i = 0; i < N; i++) x[i] += alfa * q[i];
        for (i = 0; i < N; i++) R1[i] = R[i] - alfa * Aq[i];
        double betta1 = 0;
        for (i = 0; i < N; i++) betta1 += R1[i] * R1[i];
        double betta = betta1 / alfa1;
        double u1 = 0, u2 = 0;
        for (i = 0; i < N; i++) {
            q[i] = R1[i] + betta * q[i];
            R[i] = R1[i];
            u1 += R1[i] * R1[i];
            u2 += b[i] * b[i]; }
        nrm = sqrt(u1) / sqrt(u2);
    }
}


/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void timovkinayu::lab8()
{
double E = 1.0e-6; double errc = 0;
    double c = 0, s = 0, alpha = 0; int MI = 0, MJ = 1;
    double** C; C = new double* [N];
    for (int i = 0; i < N; i++)
    {
        C[i] = new double[N];
        for (int j = i + 1; j < N; j++)
        {
            if (i != j) errc += A[i][j] * A[i][j];
            C[i][j] = 0; }
    }
    while (sqrt(errc) > E)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = i + 1; j < N; j++)
            {
                if (abs(A[i][j]) >= abs(A[MI][MJ])){ MI = i; MJ = j; }
            }
        }
        if (A[MI][MI] == A[MJ][MJ]) alpha = 3.14151 / 4;
        else alpha = 0.5 * atan((2 * A[MI][MJ]) / (A[MJ][MJ] - A[MI][MI]));
        c = cos(alpha); s = sin(alpha);
        C[MI][MI] = c * c * A[MI][MI] - 2 * s * c * A[MI][MJ] + s * s * A[MJ][MJ];
        C[MJ][MJ] = s * s * A[MI][MI] + 2 * s * c * A[MI][MJ] + c * c * A[MJ][MJ];
        C[MI][MJ] = (c * c - s * s) * A[MI][MJ] + s * c * (A[MI][MI] - A[MJ][MJ]);
        C[MJ][MI] = C[MI][MJ];
        for (int k = 0; k < N; k++)
        {
            if (k != MI && k != MJ)
            {
                C[MI][k] = c * A[MI][k] - s * c * A[MJ][k];
                C[k][MI] = C[MI][k];
                C[MJ][k] = s * A[MI][k] + c * A[MJ][k];
                C[k][MJ] = C[MJ][k]; }
            for (int l = 0; l < N; l++) { if (k != MI && k != MJ && l != MI && l != MJ)C[k][l] = A[k][l];}
        }
        errc = 0;
        for (int i = 0; i < N; i++)
        {
            for (int j = i + 1; j < N; j++) if (i != j) { errc += C[i][j] * C[i][j]; }
        }
        for (int i = 0; i < N; i++)
        { for (int j = 0; j < N; j++) A[i][j] = C[i][j]; }
        for (int i = 0; i < N; i++) x[i] = A[i][i];
    }
}


/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void timovkinayu::lab9()
{
double eps = 1.0e-3;
    double *ukone, *uktwo;
    double lambdaKONE = 0, lambdaKTWO = 1;
    int n = N;
    uktwo = new double[n];
    ukone = new double[n];
    for (int i = 0; i < n; i++) uktwo[i] = 1;
    while (fabs(lambdaKTWO - lambdaKONE) > eps) {
        lambdaKONE = lambdaKTWO;
        for (int i = 0; i < n; i++) {
            double s = 0;
            for (int j = 0; j < n; j++) s += A[i][j] * uktwo[j];
            ukone[i] = s;
        }
        for (int i = 0; i < n; i++) {
            if (ukone[0] != 0 && uktwo[0] != 0) lambdaKTWO = ukone[i] / uktwo[i];
        }
        for (int i = 0; i < n; i++) uktwo[i] = ukone[i];
    }
    cout << "lambda: " << lambdaKTWO << endl;
}


std::string timovkinayu::get_name()
{
  return "A.YU. Timovkin";
}

#include "melyakinev.h"

void melyakinev::lab1() // введение в дисциплину
{
  cout << "Hello World!" << endl;
}

void melyakinev::lab2() // метод гаусса с выбором главного элемента
{
    // приведение к треугольному виду
    for (int i=0; i<N; i++) {
        // выбор главного элемента
        int m = i;
        for (int j=i+1; j<N; j++)
            if (abs(A[m][i])<abs(A[j][i]))
                m = j;
        swap(A[i],A[m]);
        swap(b[i],b[m]);
        // нормировка
        for (int j=i+1; j<N; j++)
            A[i][j] /= A[i][i];
        b[i] /= A[i][i];
        A[i][i] = 1;
        // исключение столбца
        for (int j=i+1; j<N; j++) {
            for (int k=i+1; k<N; k++)
                A[j][k] -= A[i][k]*A[j][i];
            b[j] -= b[i]*A[j][i];
            A[j][i] = 0;
        }
    }
    // обратный ход
    for (int i=N-1; i>=0; i--)
        for (int j=i-1; j>=0; j--) {
            b[j] -= b[i]*A[j][i];
            A[j][i] = 0;
        }
    // ответ
    for (int i=0; i<N; i++)
        x[i] = b[i];
}

void melyakinev::lab3() // метод прогонки
{
    double* alpha = new double[N];
    double* beta = new double[N];
    // вычисление коэфициентов
    alpha[0] = A[0][1]/A[0][0];
    beta[0] = b[0]/A[0][0];
    for (int i=1; i<N; i++) {
        alpha[i] = A[i][i+1]/(A[i][i]-A[i][i-1]*alpha[i-1]);
        beta[i] = (b[i]-A[i][i-1]*beta[i-1])/(A[i][i]-A[i][i-1]*alpha[i-1]);
    }
    // ответ
    x[N-1] = beta[N-1];
    for (int i=N-2; i>=0; i--)
        x[i] = beta[i]-alpha[i]*x[i+1];
    delete[] alpha;
    delete[] beta;
}

void melyakinev::lab4() // метод холецкого
{
    double** S = new double*[N];
    for (int i = 0; i < N; i++)
        S[i] = new double[N];
    double* D = new double[N];
    for (int i = 0; i < N; i++) {
        for (int l = 0; l < i; l++)
            A[i][i] -= S[l][i] * S[l][i] * D[l];
        D[i] = A[i][i] >= 0 ? 1 : -1;
        S[i][i] = sqrt(D[i] * A[i][i]);
        for (int j = i + 1; j < N; j++) {
            for (int l = 0; l < j; l++)
                A[i][j] -= S[l][i] * D[l] * S[l][j];
            S[i][j] = A[i][j] / (S[i][i] * D[i]);
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++)
            b[i] -= S[j][i] * b[j] * D[j];
        b[i] /= S[i][i] * D[i];
    }
    for (int i = N - 1; i >= 0; i--) {
        for (int j = i + 1; j < N; j++)
            b[i] -= S[i][j] * x[j];
        x[i] = b[i] / S[i][i];
    }
    for (int i = 0; i < N; i++)
        delete[] S[i];
    delete[] S;
    delete[] D;
}

void melyakinev::lab5() // метод якоби
{
    double EPS = 1e-15;
    double* newX = new double[N];
    bool f;
    do {
        f = false;
        for (int i = 0; i < N; i++) {
            newX[i] = b[i];
            for (int j = 0; j < i; j++)
                newX[i] -= A[i][j] * x[j];
            for (int j = i + 1; j < N; j++)
                newX[i] -= A[i][j] * x[j];
            newX[i] /= A[i][i];
            if (abs((newX[i] - x[i])) > EPS)
                f = true;
        }
        for (int i = 0; i < N; i++)
            x[i] = newX[i];
    } while (f);
    delete[] newX;
}

void melyakinev::lab6() // метод минимальных невязок
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
            for (int j = 0; j < n; j++)
                r[i] += A[i][j] * x[j];
        }
        double lower = 0, upper = 0;
        for (int i = 0; i < n; i++) {
            double temp = 0;
            for (int j = 0; j < n; j++)
                temp += A[i][j] * r[j];
            lower += temp * temp;
            upper += temp * r[i];
        }
        double t = upper / lower;
        for (int i = 0; i < n; i++) {
            if (z < abs(t * r[i]))
                z = abs(t * r[i]);
            x[i] -= t * r[i];
        }
    }
}

void melyakinev::lab7() // метод сопряженных градиентов
{

}

void melyakinev::lab8() // метод вращения для нахождения собственных значений матрицы
{

}

void melyakinev::lab9() // нахождение наибольшего по модулю собственного значения матрицы
{

}

std::string melyakinev::get_name()
{
  return "E.V. Melyakin";
}

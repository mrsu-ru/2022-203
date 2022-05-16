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
    double* m = new double[N];
    double* n = new double[N];
    // вычисление коэфициентов
    m[0] = A[0][1]/A[0][0];
    n[0] = b[0]/A[0][0];
    for (int i=1; i<N; i++) {
        m[i] = A[i][i+1]/(A[i][i]-A[i][i-1]*m[i-1]);
        n[i] = (b[i]-A[i][i-1]*n[i-1])/(A[i][i]-A[i][i-1]*m[i-1]);
    }
    // ответ
    x[N-1] = n[N-1];
    for (int i=N-2; i>=0; i--)
        x[i] = n[i]-m[i]*x[i+1];
    delete[] m;
    delete[] n;
}

void melyakinev::lab4() // метод холецкого
{
    
}

void melyakinev::lab5() // метод якоби или зейделя
{

}

void melyakinev::lab6() // метод минимальных невязок
{

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

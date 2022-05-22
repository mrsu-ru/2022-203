#pragma once
#include "lab.h"

void method_Gauss(double **, double *, int);
void MatrixTrans(double **, double **, int);
void init_matrix(double *, int);
void init_matrix(double **, int);
void create_dynamic_array(double **, int);
void delete_dynamic_array(double **, int);
double *MatrixMultOnVect(double **, double *, int);
double **MatrixMultOnMatrix(double **, double **, int);
double scalar(double *, double *, int);

class melkonyanma : public lab
{

  /**
   * Введение
   */
  virtual void lab1();

  /**
   * Метод Гаусса с выбором главного элемента
   */
  virtual void lab2();

  /**
   * Метод прогонки
   */
  virtual void lab3();

  /**
   * Метод квадратного корня (метод Холецкого)
   */
  virtual void lab4();

  /**
   * Метод Якоби, Зейделя
   */
  virtual void lab5();

  /**
   * Метод минимальных невязок
   */
  virtual void lab6();

  /**
   * Метод сопряженных градиентов
   */
  virtual void lab7();

  /**
   * Метод вращения для нахождения собственных значений матрицы
   */
  virtual void lab8();

  /**
   * Нахождение наибольшего по модулю собственного значения матрицы
   */
  virtual void lab9();

  void method_Gauss(double **A, double *b, int N)
  {
    for (int row = 0, column = 0; column < N; ++row, ++column)
    {
      int max = row;

      for (int i = row + 1; i < N; ++i)
      {
        if (fabs(A[i][column]) > fabs(A[max][column]))
        {
          max = i;
        }
      }

      if (max != row)
      {
        swap(A[row], A[max]);
        swap(b[row], b[max]);
      }

      double mainElem = A[row][column];

      for (int i = 0; i < N; ++i)
      {
        if (A[row][i] == 0)
        {
          continue;
        }
        A[row][i] /= mainElem;
      }

      b[row] /= mainElem;

      for (int i = row + 1; i < N; ++i)
      {
        double a = A[i][column];
        for (int j = 0; j < N; ++j)
        {
          A[i][j] -= A[row][j] * a;
        }
        b[i] -= b[row] * a;
      }
    }

    for (int i = N - 1; i > 0; --i)
    {
      double a = A[i][i];
      for (int k = i - 1; k >= 0; --k)
      {
        b[k] -= b[i] * A[k][i];
        A[k][i] -= a * A[k][i];
      }
    }
  }

  void MatrixTrans(double **matrix, double **tMatrix, int N)
  {
    for (int i = 0; i < N; ++i)
    {
      for (int j = 0; j < N; ++j)
      {
        tMatrix[i][j] = matrix[j][i];
      }
    }
  }

  double *MatrixMultOnVect(double **A, double *V, int N)
  {
    double *res = new double[N];
    for (int i = 0; i < N; ++i)
    {
      res[i] = 0;
      for (int j = 0; j < N; ++j)
      {
        res[i] += A[i][j] * V[j];
      }
    }
    return res;
  }

  double **MatrixMultOnMatrix(double **first_m, double **second_m, int N)
  {
    double **res = new double *[N];
    create_dynamic_array(res, N);
    for (int i = 0; i < N; ++i)
    {
      for (int j = 0; j < N; ++j)
      {
        res[i][j] = 0;
        for (int k = 0; k < N; ++k)
        {
          res[i][j] += first_m[i][k] * second_m[k][j];
        }
      }
    }
    return res;
  }

  double scalar(double *A, double *B, int N)
  {
    double res = 0;
    for (int i = 0; i < N; ++i)
    {
      res += A[i] * B[i];
    }
    return res;
  }

  void init_matrix(double *matrix, int N)
  {
    for (int i = 0; i < N; ++i)
    {
      matrix[i] = 0.;
    }
  }

  void init_matrix(double **matrix, int N)
  {
    for (int i = 0; i < N; ++i)
    {
      for (int j = 0; j < N; ++j)
      {
        matrix[i][j] = 0.;
      }
    }
  }

  void create_dynamic_array(double **array, int N)
  {
    for (int i = 0; i < N; ++i)
    {
      array[i] = new double[N];
    }
  }

  void delete_dynamic_array(double **array, int N)
  {
    for (int i = 0; i < N; ++i)
    {
      delete[] array[i];
    }
  }

  virtual std::string get_name();
};

#pragma once
#include "lab.h"


class nikishkinev : public lab
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


  virtual std::string get_name();

};

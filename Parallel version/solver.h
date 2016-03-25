#ifndef _SOLVER_H_
#define _SOLVER_H_
#include <iostream>
#include <cmath>
#include <cassert>
#include <math.h>

class solver {
  public:
    solver();
    ~solver();
    virtual double rootFinder(double (*f)(double, double*, int, int, double, int, double),
                              double lb, double ub, double* s, int n, int k, double T, int r, double c);
    void set_tol(double t);
    void set_max_iter(int n);
  protected:
    double tol = 1e-6;
    int max_iter = 150;
};

class bisection: public solver {
  public:
    double rootFinder(double (*f)(double, double*, int, int, double, int, double),
                      double lb, double ub, double* s, int n, int k, double T, int r, double c);
};

class newton: public solver {
  public:
    double rootFinder(double (*f)(double, double*, int, int, double, int, double),
                      double lb, double ub, double* s, int n, int k, double T, int r, double c);
};

class brent: public solver {
  public:
    double rootFinder(double (*f)(double, double*, int, int, double, int, double),
                      double lb, double ub, double* s, int n, int k, double T, int r, double c);
};
#endif



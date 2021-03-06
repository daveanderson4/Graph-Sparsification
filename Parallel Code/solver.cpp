#include "solver.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>

solver::solver() {};

solver::~solver() {};

double solver::rootFinder(double (*f)(double, std::vector<double>&, int, int, double, int, double),
                          double lb, double ub, std::vector<double>& s, int n, int k, double T, int r,
			  double c) { return 0; };

void solver::set_tol(double t) { tol = t; };

void solver::set_max_iter(int n) { max_iter = n; };

double bisection::rootFinder(double (*f)(double, std::vector<double>&, int, int, double, int, double),
                             double lb, double ub, std::vector<double>& s, int n, int k, double T, int r,
			     double c) {
  // find root with bisection method
  //   - stable, but slow

  set_max_iter(1000);
  int iter = 0;
  double mid = 0.5*(lb+ub);
  while (iter<max_iter && std::abs((*f)(mid, s, n, k, T, r, c)) > tol) {
    (*f)(lb, s, n, k, T, r, c) * (*f)(mid, s, n, k, T, r, c) < 0 ? ub = mid : lb = mid;
    mid = 0.5*(lb+ub);
    ++iter;
  }
  assert(std::abs((*f)(mid, s, n, k, T, r, c)) <= tol); // method reached max_iter without converging
  return mid;
};

double newton::rootFinder(double (*f)(double, std::vector<double>&, int, int, double, int, double),
                          double lb, double ub, std::vector<double>& s, int n, int k, double T, int r,
			  double c) {
  // find root with newton's method
  //   - fast, but unstable
  set_max_iter(30);
  std::cout << "newton not implemented" << std::endl;
  return 0.0;
};

double brent::rootFinder(double (*f)(double, std::vector<double>&, int, int, double, int, double),
                         double lb, double ub, std::vector<double>& s, int n, int k, double T, int r,
			 double c) {
  // find root with brent zero method
  //   - stability of bisection, with speed of newton's
  set_max_iter(30);
  std::cout << "brent not implemented" << std::endl;
  return 0.0;
};




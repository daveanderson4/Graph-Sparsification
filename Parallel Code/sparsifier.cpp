#include "sparsifier.h"
#include "solver.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <float.h>
#include <vector>
#include <cassert>
#include <memory>

// link to LAPACK functions
extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

    // copy a matrix
    void dlacpy_(char* UPLO, int* M, int* N, double* A, int* LDA, double* B, int* LDB);

    // matrix-vector multiply
    void dgemv_(char* TRANS, int* M, int* N, double* ALPHA, double* A, int* LDA, 
                double* X, int* INCX, double* BETA, double* Y, int* INCY);

    // find eigenvalues
    void dgeev_(char* JOBVL, char* JOBVR, int* N, double* A, int* LDA, double* WR, 
                double* WI, double* VL, int* LDVL, double* VR, int* LDVR, double* WORK,
                int* LWORK, int* INFO);

    // symmetric eigenvalue problem (replaces dgeev)
    void dsyev_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA, double* W, 
                double* WORK, int* LWORK, int* INFO);

    // symmetric eigenvalue problem (replaces dgeev)
    void dsyevd_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA, double* W,
                double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);
}

sparsifier::sparsifier(int n, int k, int l) : n(n), k(k), l(l) {
  assert(n>=l && l>k && k>0);
  T = findT(n,k,l);
};

sparsifier::~sparsifier() {};

void sparsifier::columnSelect(int* pi, double* U) {
  std::cout << "columnSelect is defunct" << std::endl;
};

int sparsifier::getN() { return n; };

double sparsifier::getT() { return T; };

double sparsifier::findT(int n, int k, int l) {
  // finds appropriate progress parameter T
  double nd = (double)n;
  double kd = (double)k;
  double ld = (double)l;
  double T_hat = (kd*(nd+(ld+1)/2-kd)+sqrt(kd*ld*(nd-(ld-1)/2)
                    *(nd+(ld+1)/2-kd)))/(ld-kd);
  double T = T_hat*(1+(1-kd/T_hat)*ld/(nd-(ld-1)/2-kd+T_hat)-kd/T_hat);
  return T;
};

double sparsifier::root1(double x, std::vector<double>& s, int n, int k, double T, int r, double c) {
  // first root function
  double root = 0;
  for (auto si : s) {
    root += 1 / (si - x);
  }
  return (root - T);
};

double sparsifier::root2helper(std::vector<double>& s, int n, int k, double lambda, int r) {
  double c = 0;
  for (auto si : s) {
    c += (1 - si) / (si - lambda);
  }
  c += (double)n - (double)r;
  return c;
};

double sparsifier::root2(double lh, std::vector<double>& s, int n, int k, double lambda, int r, double c) {
  // second root function
  double root = 0;
  root += (lh - lambda) * c;
  double n1 = 0, n2 = 0;
  for (auto si : s) {
    n1 += (1 - si) / ((si - lambda) * (si - lh));
    n2 += 1 / ((si - lambda) * (si - lh));
  }
  root -= (n1 / n2);
  return root;
};

void sparsifier::inverse(std::shared_ptr<double> A, int N) {
  // for the inverse of a N-by-N matrix
  int INFO, LWORK = N * N;
  std::shared_ptr<int>    IPIV (new int[N + 1]);
  std::shared_ptr<double> WORK (new double[LWORK]);
  dgetrf_(&N, &N, A.get(), &N, IPIV.get(), &INFO);
  dgetri_(&N, A.get(), &N, IPIV.get(), WORK.get(), &LWORK, &INFO);
};












#include "sparsifier.h"
#include "solver.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <float.h>
#include <cassert>

#ifdef __APPLE__
  #include <Accelerate/Accelerate.h>
#else
  // link to LAPACK functions
  extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* LDA, int* IPIV, double* WORK, int* LWORK, int* INFO);

    // copy a matrix
    void dlacpy_(char* UPLO, int* M, int* N, double* A, int* LDA, double* B, int* LDB);

    // matrix-vector multiply
    void dgemv_(char* TRANS, int* M, int* N, double* ALPHA, double* A, int* LDA, 
                double* X, int* INCX, double* BETA, double* Y, int* INCY);

    // symmetric eigenvalue problem
    void dsyev_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA, double* W, 
                  double* WORK, int* LWORK, int* INFO);
  }
#endif

using namespace std;

void sparsifiers_dgemv(double* Ainv, int k, int n, double* space, double* space2, int INCX, int j) {
  #ifdef __APPLE__
    cblas_dgemv(CblasColMajor, CblasNoTrans, k, k, 1., Ainv, k, &space2[j], INCX, 0., space, 1);
  #else
    char TRANS = 'N';
    double ONE = 1.0, ZERO = 0.0;
    int ONEINT = 1;
    dgemv_(&TRANS, &k, &k, &ONE, Ainv, &k, &space2[j], &INCX, &ZERO, space, &ONEINT);
  #endif
};

void sparsifiers_copy(double* A, double* Acopy, int k) {
  #ifdef __APPLE__
    cblas_dcopy(k*k, A, 1, Acopy, 1);
  #else
    char UPLO = 'N';
    dlacpy_(&UPLO, &k, &k, A, &k, Acopy, &k);
  #endif
};

sparsifier::sparsifier(int n, int k, int l) : n(n), k(k), l(l) {
  assert(n>=l && l>k && k>0);
  T = findT(n,k,l);
};

sparsifier::~sparsifier() {};

void sparsifier::columnSelect(int* pi, double* U) {
  // pi should be all zeros, length n
  double lambda     = -DBL_MAX; // smallest double
  double lambda_hat = -DBL_MAX;
  int    mindx = -1;
  int INFO;
  double* A = new double[k*k];
  double* space = new double[k];
  bisection b1;
  solver* s1 = &b1;
  for (int i=0; i<k*k; ++i) {A[i] = 0;}
  double slambda;
  double slambdah;
  slambda=lambda;
  slambdah = lambda_hat;
 
  // iterate
  for (int r=0; r<l; ++r) {
    lambda = slambda;
    lambda_hat = slambdah;
    char    NO = 'N';
    double* Acopy = new double[k*k];
    double* s     = new double[k];
    double* WI    = new double[k];

    sparsifiers_copy(A, Acopy, k);

    char UP = 'U'; // both are stored
    int LWORK = 3*k-1;  // optimal size
    double* WORK = new double[LWORK];

    dsyev_(&NO, &UP, &k, Acopy, &k, s, WORK, &LWORK, &INFO);
 
    delete Acopy;
    delete WI;
    delete WORK;

    double mins = max(0., s[0]);
    double r1_lb = max(lambda_hat, mins-1);
    double r1_ub = mins-1e-12;

    assert(root1(r1_lb, s, n, k, T, 0, 0.)*root1(r1_ub, s, n, k, T, 0, 0.) < 0);
    double (*root)(double, double*, int, int, double, int, double) = &sparsifier::root1;
    double lambda = s1->rootFinder(root, r1_lb, r1_ub, s, n, k, T, 0, 0.);

    // step 2: compute lambda_hat
    double r2_lb = lambda;
    double r2_ub = mins - 1e-12;
    double c = root2helper(s, n, k, lambda, r);

    assert(root2(r2_lb, s, n, k, lambda, r, c)*root2(r2_ub, s, n, k, lambda, r, c) < 0);
    double (*secroot)(double, double*, int, int, double, int, double) = &sparsifier::root2;
    double lambda_hat = s1->rootFinder(secroot, r2_lb, r2_ub, s, n, k, lambda, r, c);
    slambda = lambda;
    slambdah = lambda_hat;

    // step 3: look for index
    double trinv=0, trhinv=0;
    for (int i=0; i<k; i++) {
      trinv  += 1/(s[i] - lambda);
      trhinv += 1/(s[i] - lambda_hat);
    }
    delete s;
    double* Ainv = new double[k*k];

    sparsifiers_copy(A, Ainv, k);

    for (int ii=0; ii<k; ++ii) { Ainv[(k+1)*ii] -= lambda_hat; }

    inverse(Ainv, k);

    for (int i=0; i<n; ++i) {
      int j = ((i+mindx+1) % n); // current search index
      if (pi[j] == 0) {          // not already selected
        for (int jj=0; jj<k; ++jj) {
          space[jj] = 0;
          for (int m=0; m<k; ++m) {
            space[jj] += U[n*m+j] * Ainv[jj+k*m];
          }
        }
        double d = 1.0;
        for (int jj=0; jj<k; ++jj) {
          d += space[jj]*U[n*jj+j];
        }
        for (int ll=0; ll<k; ++ll) { space[ll] = U[n*ll+j]; }

        sparsifiers_dgemv(Ainv, k, n, space, U, n, j);

        double* space2 = new double[k];
        for (int pp=0; pp<k; ++pp) { space2[pp] = space[pp]; }

        sparsifiers_dgemv(Ainv, k, n, space, space2, 1, 0);

        delete space2;
        double adj = 0;
        for (int ii=0; ii<k; ++ii) { adj += U[j+n*ii]*space[ii]; }
        adj /= d;
          
        // test column j
        if (trhinv-adj <= trinv) {
          for (int ii=0; ii<k; ++ii) {
            for (int y=0; y<k; ++y) {
              A[k*y+ii] += U[j+n*ii]*U[j+n*y];
            }
          }
          pi[j] = r + 1;
          mindx = j; // save where loop ended
          break;
        }
      }
    }
    delete Ainv;
  }
  delete A;
  delete space;
};

int sparsifier::getN() {return n;};

double sparsifier::getT() {return T;};

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

double sparsifier::root1(double x, double* s, int n, int k, double T, int r, double c) {
  // first root function
  double root = 0;
  for (int i=0; i<k; ++i) {
    root += 1/(s[i]-x);
  }
  return (root - T);
};

double sparsifier::root2helper(double* s, int n, int k, double lambda, int r) {
  double c = 0;
  for (int i=0; i<k; ++i) {
    c += (1-s[i])/(s[i]-lambda);
  }
  c += (double)n-(double)r;
  return c;
};

double sparsifier::root2(double lh, double* s, int n, int k, double lambda, int r, double c) {
  // second root function
  double root = 0;
  root += (lh-lambda)*c;
  double n1 = 0, n2 = 0;
  for (int i=0; i<k; ++i) {
    n1 += (1-s[i])/((s[i]-lambda)*(s[i]-lh));
    n2 += 1/((s[i]-lambda)*(s[i]-lh));
  }
  root -= (n1/n2);
  return root;
};

void sparsifier::inverse(double* A, int N) { 
  // for the inverse of a N-by-N matrix
  int *IPIV = new int[N];
  int LWORK = N*N;
  double *WORK = new double[LWORK];
  int INFO;
  dgetrf_(&N,&N,A,&N,IPIV,&INFO);
  dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
  delete IPIV;
  delete WORK;
};








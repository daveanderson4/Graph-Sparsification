#include "sparsifier.h"
#include "solver.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <float.h>
#include <cassert>

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

using namespace std;

sparsifier::sparsifier(int n, int k, int l) : n(n), k(k), l(l) {
  assert(n>=l && l>k && k>0);
  T = findT(n,k,l);
};

sparsifier::~sparsifier() {};

void sparsifier::columnSelect(int* pi, double* U) {
  cout << "columnSelect is defunct" << endl;
  /*
  // pi should be all zeros, length n
  double lambda     = -1e308;//-DBL_MAX; // smallest double
  double lambda_hat = -1e308;//-DBL_MAX;
  int    mindx = -1;
  double lbond3 = -(double)n/T;
  char UPLO = 'N';
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
    double* Acopy = new double[k*k]; // eigenvalue routine destroys matrix
    double* s     = new double[k];
    double* WI    = new double[k];
    int     LWORK = 20*k; // LWORK >= (NB+2)*k optimally
    double* WORK  = new double[LWORK];
    double *VL, *VR;
    int     ONEINT = 1;
    dlacpy_(&UPLO, &k, &k, A, &k, Acopy, &k);

    char UP = 'U'; // both parts are stored in this implementation
    int LIWORK = k;
    int* IWORK = new int[LIWORK];
    dsyevd_(&NO, &UP, &k, Acopy, &k, s, WORK, &LWORK, IWORK, &LIWORK, &INFO);
    delete IWORK;

    delete Acopy;
    delete WI;
    delete WORK;

    double mins = max(s[0], 0.); // dsyev reports eigenvalues in ascending order

    double r1_lb = max(lambda_hat, mins-1);
    double r1_ub = mins-1e-12;

    assert(root1(r1_lb, s, n, k, T, 0)*root1(r1_ub, s, n, k, T, 0) < 0);
    double (*root)(double, double*, int, int, double, int) = &sparsifier::root1;
    double lambda = s1->rootFinder(root, r1_lb, r1_ub, s, n, k, T, 0);

    // step 2: compute lambda_hat
    double r2_lb = lambda;
    double r2_ub = mins - 1e-12;
    assert(root2(r2_lb, s, n, k, lambda, r)*root2(r2_ub, s, n, k, lambda, r) < 0);
    double (*secroot)(double, double*, int, int, double, int) = &sparsifier::root2;
    double lambda_hat = s1->rootFinder(secroot, r2_lb, r2_ub, s, n, k, lambda, r);
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
    dlacpy_(&UPLO, &k, &k, A, &k, Ainv, &k);
    for (int ii=0; ii<k; ++ii) { Ainv[(k+1)*ii] -= lambda_hat; }
    inverse(Ainv, k);
    for (int i=0; i<n; ++i) {
      int j = ((i+mindx+1) % n); // current search index
      if (pi[j] == 0) { // not already selected
        for (int jj=0; jj<k; ++jj) {
          space[jj] = 0;
          for (int m=0; m<k; ++m) {
            double entry = Ainv[jj+k*m];
            //double diag = i==jj ? -lambda_hat : 0;
            space[jj] += U[n*m+j]*entry;//(entry+diag);
          }
        }
        double d = 1.0;
        for (int jj=0; jj<k; ++jj) {
          d += space[jj]*U[n*jj+j];
        }
        for (int ll=0; ll<k; ++ll) { space[ll] = U[n*ll+j]; }
        char TRANS = 'N';
        double ONE = 1.0, ZERO = 0.0;
        int INCX = 1;
        dgemv_(&TRANS, &k, &k, &ONE, Ainv, &k, &U[j], &n,    &ZERO, space, &INCX);
        double* space2 = new double[k];
        for (int pp=0; pp<k; ++pp) { space2[pp] = space[pp]; }
        dgemv_(&TRANS, &k, &k, &ONE, Ainv, &k, space2, &INCX, &ZERO, space, &INCX);
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
    Ainv = NULL;
    delete Ainv;
  }
  A = NULL;
  space = NULL;
  delete A;
  delete space;
  */
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
  int *IPIV = new int[N+1];
  int LWORK = N*N;
  double *WORK = new double[LWORK];
  int INFO;

  dgetrf_(&N,&N,A,&N,IPIV,&INFO);
  dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

  IPIV = NULL;
  WORK = NULL;
  delete IPIV;
  delete WORK;
};












#include "sparsifier.h"
#include "solver.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <float.h>
#include <stdlib.h>
#include <cassert>
#include <mpi.h>
#include <sys/time.h>

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

double read_timer() {
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

void convertPi(int* pi_in, int* pi_out, int n, int k) {
  for (int i=0; i<n; ++i) {
    if (pi_in[i] != 0) {
      pi_out[pi_in[i]-1] = i+1;
    }
  }
}

int main(int argc, char* argv[]) {

  // initialize MPI
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // the initial work will be carried out on the root processor
  bool root = (rank == 0);

  // declare parameters
  int n, k, l;
  double lambda, slambda, lambda_hat, slambdah;
  bool VERBOSE, WRITE;
  char UPLO;
  int INFO;
  double *space, *A, *U, *Ainv;
  int *pi, *pi2;
  double T;
  int mindx;
  double trinv=0, trhinv=0;
  solver *s1;
  bisection b1;
  bool record = true;
  double seconds;

  // parameters
  n = 73;
  k = 35;
  l = 2*k;
  if (argc > 3) {
    n = atoi(argv[1]);
    k = atoi(argv[2]);
    l = atoi(argv[3]);
  }
  sparsifier S(n, k, l);
  T = S.getT();

  if (root) {
    VERBOSE = false;
    WRITE   = true;
    
    // read in data
    ifstream fin;
    fin.open ("U.txt");
    assert(fin);
    U = new double[n*k];
    for (int i=0; i<n; ++i) {
      for (int j=0; j<k; ++j) {
        fin >> U[i+n*j];
      }
    }
    fin.close();

    // begin timer here
    if (root) seconds = read_timer();

    pi = new int[n];
    for (int i=0; i<n; ++i) { pi[i] = 0; }
  
    lambda     = -1e308;//-DBL_MAX; // smallest double
    lambda_hat = -1e308;//-DBL_MAX;
    mindx = -1;

    double lbond3 = -(double)n/T;
    UPLO = 'N';
    A = new double[k*k];
    space = new double[k];

    s1 = &b1;

    for (int i=0; i<k*k; ++i) {A[i] = 0;}
    slambda=lambda;
    slambdah = lambda_hat;
  } // ends if(root)

  // all processors need (part) of U
  if (!root) U = new double[n*k];
  MPI_Bcast(U,  n*k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  // iterate
  for (int r=0; r<l; ++r) {

    // perform easy calculations on root only
    if (root) {

      // superfluous send on previous iteration causes a blank pass through
      if (record){

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

        assert(S.root1(r1_lb, s, n, k, T, 0, 0.)*S.root1(r1_ub, s, n, k, T, 0, 0.) < 0);
        double (*firstroot)(double, double*, int, int, double, int, double) = &sparsifier::root1;
        double lambda = s1->rootFinder(firstroot, r1_lb, r1_ub, s, n, k, T, 0, 0.);

        // step 2: compute lambda_hat
        double r2_lb = lambda;
        double r2_ub = mins - 1e-12;

        double c = S.root2helper(s, n, k, lambda, r);

        assert(S.root2(r2_lb, s, n, k, lambda, r, c)*S.root2(r2_ub, s, n, k, lambda, r, c) < 0);
        double (*secroot)(double, double*, int, int, double, int, double) = &sparsifier::root2;
        double lambda_hat = s1->rootFinder(secroot, r2_lb, r2_ub, s, n, k, lambda, r, c);
        slambda = lambda;
        slambdah = lambda_hat;

        // step 3: look for index
        // search is performed in parallel
        trinv=0, trhinv=0;
        for (int i=0; i<k; i++) {
          trinv  += 1/(s[i] - lambda);
          trhinv += 1/(s[i] - lambda_hat);
        }
        delete s;

        Ainv = new double[k*k];
        dlacpy_(&UPLO, &k, &k, A, &k, Ainv, &k);
        for (int ii=0; ii<k; ++ii) { Ainv[(k+1)*ii] -= lambda_hat; }
        S.inverse(Ainv, k);

      } // end if(record)
    } // end if(root)

    // other cores need to allocate space
    if (!root) {
      Ainv = new double[k*k];
      space = new double[k];
      pi = new int[n];
    }

    // send variables out to other cores
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(Ainv,  k*k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&trinv,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&trhinv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(pi,      n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // begin unblocked receive
    int fcol;
    MPI_Status status;
    MPI_Request request;
    MPI_Irecv(&fcol, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &request);

    // test if done already
    int done = 0;
    MPI_Test(&request, &done, &status);
    
    // find bound on search indices for processor
    int start  = n/size*rank;     // inclusive
    int finish = n/size*(rank+1); // exclusive
    if (rank == size-1) finish = n;

    MPI_Barrier(MPI_COMM_WORLD);

    // if this iteration did not find a new index then don't record anything!
    record = !done;
    // if this iteration did nothing then the work remains to be computed
    if (!record) --r;

    for (int j=start; j<finish && !done; ++j) {
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
          for (int prc=0; prc<size; ++prc) {
            fcol = j;
            MPI_Send(&fcol, 1, MPI_INT, prc, 1, MPI_COMM_WORLD);
          }
        }
      }

      // test if found
      MPI_Test(&request, &done, &status);

    } // end of search loop

    if (!root) {
      Ainv = NULL;
      //delete Ainv;
      space = NULL;
      //delete space;
      pi = NULL;
      //delete pi;
    }

    MPI_Wait(&request, &status);
    MPI_Barrier(MPI_COMM_WORLD);

    // update A
    if (root) {

      // only record if this iteration actually found a new entry
      if (record) {
        for (int ii=0; ii<k; ++ii) {
          for (int y=0; y<k; ++y) {
            A[k*y+ii] += U[fcol+n*ii]*U[fcol+n*y];
          }
        }
        pi[fcol] = r+1;
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

  } // end of iteration for loop
  if (root) {
    A = NULL;
    space = NULL;
    delete A;
    delete space;
  }
  // column selection is complete

  if (root) {
    // record selection
    if (VERBOSE) {
      cout << "pi: " << endl;
      for (int i=0; i<n; ++i) { cout << pi[i] << endl; }
      cout << endl;
    }
    pi2 = new int[l];
    convertPi(pi, pi2, n, l);
    if (VERBOSE) {
      cout << "pi2: " << endl;
      for (int i=0; i<l; ++i) { cout << pi2[i] << endl; }
      cout << endl;
    }
    if (WRITE) {
      ofstream pfile;
      pfile.open("p.txt");
      for (int i=0; i<l; ++i) {
        pfile << pi2[i] << "\n";
      }
      pfile.close();
    }
  }

  if (root) {
    delete U;
    delete pi;
    delete pi2;
  }

  cout <<"Processor " << rank << " has finished" << endl;

  // end timer here
  if (root) {
    seconds = read_timer() - seconds;
    cout << "Computation time: " << seconds << " seconds" << endl;
    cout << "  on " << size << " cores" << endl;
  }

  MPI_Finalize();
  return 0;
}

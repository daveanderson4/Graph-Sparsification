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
  int 		n, 	k, 	l, 	INFO, 	*pi, 	*pi2;
  double 	lambda, slambda, 	lambda_hat, 	slambdah, 
		T, 	seconds, 	trinv=0, 	trhinv=0, 
		*A, 	*Ainv, 	*U, 	*space;
  bool 		VERBOSE, 	WRITE, 	record=true;
  char 		UPLO;
  solver 	*s1;
  bisection	b1;

  // set parameters
  assert (argc > 2); // must specify problem size
  if (argc > 3) {
    n = atoi(argv[1]);
    k = atoi(argv[2]);
    if (argc > 3) {
      l = atoi(argv[3]);
    }
    else l = k+1; // default to minimum spanning tree
  }
  sparsifier S(n, k, l);
  T = S.getT();

  // find bound on search indices for processor
  int start  = n/size*rank;      	// inclusive
  int finish = n/size*(rank+1); 	// exclusive
  if (rank == size-1) finish = n;	// corner case
  int n_local = finish-start;

  // allocate space
  U = new double[n_local*k];

  double* U_dest;
  if (rank != 0) {
    U_dest = new double[n_local*k];
  }

  // sequential parts of this algorithm will be calculated on the root
  // setup information for the root here
  if (root) {
    VERBOSE = false;
    WRITE   = true;
    
    // read in data
    ifstream fin;
    fin.open ("U.txt");
    assert(fin);

    // load data for root
    for (int i=0; i<n_local; ++i) {
      for (int j=0; j<k; ++j) {
        fin >> U[i+n_local*j];
      }
    }

    // read and send data for each core
    for (int core=1; core<size; ++core) {
      int dest_start = n/size*core;
      int dest_finish = n/size*(core+1);
      if (core == size-1) dest_finish = n;
      int n_dest = dest_finish - dest_start;
      U_dest = new double[n_dest*k];
      for (int i=0; i<n_dest; ++i) {
        for (int j=0; j<k; ++j) {
          fin >> U_dest[i+n_dest*j];
        }
      }
      MPI_Send(&U_dest[0], n_dest*k, MPI_DOUBLE, core, 0, MPI_COMM_WORLD);
    }
    fin.close();

    // begin timer here
    seconds = read_timer();

    pi = new int[n];
    for (int i=0; i<n; ++i) { pi[i] = 0; }
  
    lambda     = -DBL_MAX; // smallest double
    lambda_hat = -DBL_MAX;

    double lbond3 = -(double)n/T;
    UPLO = 'N';
    A = new double[k*k];
    space = new double[k];

    s1 = &b1;

    for (int i=0; i<k*k; ++i) {A[i] = 0;}
    slambda=lambda;
    slambdah = lambda_hat;
  } // ends if(root)
  else {
    // receive data if not root
    MPI_Status recv_status;
    MPI_Recv(&U[0], n_local*k, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &recv_status);

  } // data has been sent out

  // iterate
  for (int r=0; r<l; ++r) {

    // perform easy/sequential calculations on root only
    if (root) {

      // superfluous send on previous iteration causes a blank pass through
      if (record){

        lambda = slambda;
        lambda_hat = slambdah;
        char    NO = 'N';
        int     LWORK = 20*k; // LWORK >= (NB+2)*k optimally
        int     ONEINT = 1;
        double* Acopy = new double[k*k]; // eigenvalue routine destroys matrix
        double* s     = new double[k];
        double* WI    = new double[k];
        double* WORK  = new double[LWORK];
        double *VL, *VR;
        dlacpy_(&UPLO, &k, &k, A, &k, Acopy, &k);

        char UP = 'U'; // both parts are stored in this implementation
        int LIWORK = k;
        int* IWORK = new int[LIWORK];

        dsyevd_(&NO, &UP, &k, Acopy, &k, s, WORK, &LWORK, IWORK, &LIWORK, &INFO);
        delete[] IWORK;
        delete[] Acopy;
        delete[] WI;
        delete[] WORK;

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
        delete[] s;

        Ainv = new double[k*k];
        dlacpy_(&UPLO, &k, &k, A, &k, Ainv, &k);
        for (int ii=0; ii<k; ++ii) { Ainv[(k+1)*ii] -= lambda_hat; }
        S.inverse(Ainv, k);

      } // end if(record)
    } // end if(root)

    // other cores need to allocate space
    if (!root) {
      Ainv  = new double[k*k];
      space = new double[k];
      pi    = new int[n];
    }

    // send variables out to other cores
    MPI_Bcast(Ainv,  k*k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&trinv,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&trhinv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(pi,      n, MPI_INT, 0, MPI_COMM_WORLD);

    // begin unblocked receive
    int fcol;
    MPI_Status status;
    MPI_Request request;

    // data to contain found vector, which will need to be sent from the core that found it
    // to the root
    // u_vec contains the selected column, with an additional entry indicating the location
    // of the found column
    double* u_vec;
    u_vec = new double[k+1];

    // non blocking receive
    MPI_Irecv(&u_vec[0], k+1, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &request);

    // test if done already
    int done = 0;
    MPI_Test(&request, &done, &status);

    // if this iteration did not find a new index then don't record anything!
    record = !done;
 
    // if this iteration did nothing then the work remains to be computed
    if (!record) --r;

    // search over problem subdomain for this core
    for (int j=0; j<n_local && !done; ++j) {
      int j_global = j+start;
      if (pi[j_global] == 0) { // not already selected
        for (int jj=0; jj<k; ++jj) {
          space[jj] = 0;
          for (int m=0; m<k; ++m) {
            space[jj] += U[n_local*m+j]*Ainv[jj+k*m];
          }
        } // ends for over jj

        double d = 1.0;
        for (int jj=0; jj<k; ++jj) {
          d += space[jj]*U[n_local*jj+j];
        }
        for (int ll=0; ll<k; ++ll) { space[ll] = U[n_local*ll+j]; }
        char   TRANS = 'N';
        double ONE   = 1.0, ZERO = 0.0;
        int    INCX  = 1;
        dgemv_(&TRANS, &k, &k, &ONE, Ainv, &k, &U[j], &n_local,    &ZERO, space, &INCX);
        double* space2 = new double[k];

        for (int pp=0; pp<k; ++pp) { space2[pp] = space[pp]; }
        dgemv_(&TRANS, &k, &k, &ONE, Ainv, &k, space2, &INCX, &ZERO, space, &INCX);
        delete[] space2;

        double adj = 0;
        for (int ii=0; ii<k; ++ii) { adj += U[j+n_local*ii]*space[ii]; }
        adj /= d;
 
        // test column j
        if (trhinv-adj <= trinv) {
          for (int prc=0; prc<size; ++prc) {
            // send data
            for (int i=0; i<k; ++i) u_vec[i] = U[j+n_local*i];
            u_vec[k] = j_global;
            MPI_Send(&u_vec[0], k+1, MPI_DOUBLE, prc, 1, MPI_COMM_WORLD);
          }
        }
      } // ends not already selected if statement

      // test if found
      MPI_Test(&request, &done, &status);

    } // end of search loop

    MPI_Barrier(MPI_COMM_WORLD);

    // send column (row in this formulation of U) to update A
    int find_root = fcol*size/n;

    // update A
    if (root) {

      // only record if this iteration actually found a new entry
      if (record) {
        for (int ii=0; ii<k; ++ii) {
          for (int y=0; y<k; ++y) {
            A[k*y+ii] += u_vec[ii] * u_vec[y];
          }
        }
        pi[int(u_vec[k])] = r+1;
      }
    }
  } // end of iteration for loop

  if (root) {
    delete[] A;
    delete[] space;
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
    delete[] U;
    delete[] pi;
    delete[] pi2;
  }

  cout <<"Processor " << rank << " has finished" << endl;

  // end timer here
  if (root) {
    seconds = read_timer() - seconds;
    cout << "Computation time: " << seconds << " seconds on " << size << " core(s)" << endl;
  }

  MPI_Finalize();
  return 0;
}

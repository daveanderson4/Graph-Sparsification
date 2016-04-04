#include "sparsifier.h"
#include "solver.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <float.h>
#include <stdlib.h>
#include <vector>
#include <cassert>
#include <mpi.h>
#include <memory>
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

    // symmetric eigenvalue problem (replaces dgeev)
    void dsyevd_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA, double* W,
                double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);
}

double read_timer() {
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized ) {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

void load_data(std::shared_ptr<double> U_ptr, int rank, int size, int n_local, int n, int k) {
  if (rank == 0) {
    // read in data
    std::ifstream fin;
    fin.open ("U.txt");
    assert(fin);

    // load data for root
    for (int i=0; i<n_local; ++i) {
      for (int j=0; j<k; ++j) {
        fin >> U_ptr.get()[i + n_local * j];
      }
    }

    // read and send data for each core
    for (int core=1; core<size; ++core) {
      int dest_start  = n / size * core;
      int dest_finish = n / size * (core + 1);
      if (core == size - 1) dest_finish = n;
      int n_dest = dest_finish - dest_start;
      std::shared_ptr<double> U_dest (new double [n_dest * k]);
      for (int i=0; i<n_dest; ++i) {
        for (int j=0; j<k; ++j) {
          fin >> U_dest.get()[i + n_dest * j];
        }
      }
      MPI_Send(&(U_dest.get()[0]), n_dest * k, MPI_DOUBLE, core, 0, MPI_COMM_WORLD);
    }
    fin.close();
  }
  else {
    // receive data if not root
    MPI_Status recv_status;
    MPI_Recv(&(U_ptr.get()[0]), n_local * k, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &recv_status);
  }
}

void convertPi(int* pi_in, std::shared_ptr<int> pi_out, int n, int k) {
  for (int i=0; i<n; ++i) {
    if (pi_in[i] != 0) {
      pi_out.get()[pi_in[i] - 1] = i + 1;
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

  // set parameters
  assert (argc > 2);    		// must specify problem size
  int n = atoi(argv[1]);
  int k = atoi(argv[2]);
  int l = k + 1; 			// default to minimum spanning tree
  if (argc > 3) {
    l = atoi(argv[3]);
  }
  sparsifier S(n, k, l);
  double T = S.getT();

  // find bound on search indices for processor
  int start  = n / size * rank;      	// inclusive
  int finish = n / size * (rank + 1); 	// exclusive
  if (rank == size-1) finish = n;	// corner case
  int n_local = finish - start;

  // allocate space for main data structures
  // can make changes directly to data and can pass
  // data to legacy functions
  double* A    = new double[k * k];
  double* Ainv = new double[k * k];
  double* U    = new double[n_local * k];

  // create smart pointers
  // can manage memory and send pointers
  std::shared_ptr<double> A_ptr (A);
  std::shared_ptr<double> Ainv_ptr (Ainv);
  std::shared_ptr<double> U_ptr (U);

  // other storage needed
  double *space = new double[k];
  std::vector<int> pi (n);

  // load data from file "U.txt"
  // root will read in and send data
  // other cores will wait to receive data
  load_data(U_ptr, rank, size, n_local, n, k);

  // some other variables
  double seconds;
  bool VERBOSE, WRITE, record = true;;
  solver *s1;
  bisection b1;
  double lambda, lambda_hat, trinv=0, trhinv=0;

  if (root) {
    // set some variables
    VERBOSE = false;
    WRITE   = true;

    // begin timer here
    seconds = read_timer();

    // pi must be zero to start
    for (int i=0; i<n; ++i) {
      pi[i] = 0;
    }

    // set some other parameters
    lambda     = -DBL_MAX; // smallest double
    lambda_hat = -DBL_MAX;
    s1 = &b1;

    // A must be zero to start
    for (int i=0; i<k*k; ++i) { A[i] = 0; }
  }

  // iterate
  for (int r=0; r<l; ++r) {

    // perform easy/sequential calculations on root only
    if (root) {

      // superfluous send on previous iteration causes a blank pass through
      if (record){

        // step 0: compute eigenvalues of current A before proceeding
        char UPLO = 'N';
        int  LWORK = 20*k; // LWORK >= (NB+2)*k optimally
        int  ONEINT = 1;
        std::vector<double> Acopy (k * k); // for eigenvalue routine, which destroys matrix
        dlacpy_(&UPLO, &k, &k, A, &k, &Acopy[0], &k);

        char NO = 'N';
        char UP = 'U'; // both parts are stored in this implementation
        int LIWORK = k;
        int INFO;
        std::vector<double> s (k);
        std::shared_ptr<double> WORK (new double[LWORK]);
        std::shared_ptr<int>   IWORK (new    int[LIWORK]);
        dsyevd_(&NO, &UP, &k, &Acopy[0], &k, &s[0], WORK.get(), &LWORK, IWORK.get(), &LIWORK, &INFO);
        
        // step 1: compute lambda
        // first, find interval containing lambda
        double mins  = std::max(s[0], 0.); // dsyev reports eigenvalues in ascending order
        double r1_lb = std::max(lambda_hat, mins-1);
        double r1_ub = mins-1e-12;

        // second, invoke root finder to calculate lambda
        assert(S.root1(r1_lb, s, n, k, T, 0, 0.)*S.root1(r1_ub, s, n, k, T, 0, 0.) < 0);
        double (*firstroot)(double, std::vector<double>&, int, int, double, int, double) = &sparsifier::root1;
        double lambda = s1->rootFinder(firstroot, r1_lb, r1_ub, s, n, k, T, 0, 0.);

        // step 2: compute lambda_hat
        // first, find interval containing lambda_hat and the constant c
        double r2_lb = lambda;
        double r2_ub = mins - 1e-12;
        double c= S.root2helper(s, n, k, lambda, r);

        // second, invoke root finder to calculate lambda_hat
        assert(S.root2(r2_lb, s, n, k, lambda, r, c)*S.root2(r2_ub, s, n, k, lambda, r, c) < 0);
        double (*secroot)(double, std::vector<double>&, int, int, double, int, double) = &sparsifier::root2;
        double lambda_hat = s1->rootFinder(secroot, r2_lb, r2_ub, s, n, k, lambda, r, c);

        // step 3: look for index
        // initial calculations are done on the root core
        // search is performed in parallel
        trinv=0, trhinv=0;
        for (auto s_i : s) {
          trinv  += 1/(s_i - lambda);
          trhinv += 1/(s_i - lambda_hat);
        }

        dlacpy_(&UPLO, &k, &k, A, &k, Ainv, &k);
        for (int ii=0; ii<k; ++ii) { Ainv[(k + 1) * ii] -= lambda_hat; }
        sparsifier::inverse(Ainv_ptr, k);
      } // end if(record)
    } // end if(root)

    // send variables out to other cores
    MPI_Bcast(Ainv,  k*k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&trinv,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&trhinv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&pi[0],  n, MPI_INT,    0, MPI_COMM_WORLD);

    // begin unblocked receive
    int fcol;
    MPI_Status status;
    MPI_Request request;

    // data to contain found vector, which will need to be sent from the core that found it
    // to the root
    // u_vec contains the selected column, with an additional entry indicating the location
    // of the found column
    std::vector<double> u_vec (k + 1);

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
      int j_global = j + start;
      if (pi[j_global] == 0) { // not already selected
        for (int jj=0; jj<k; ++jj) {
          space[jj] = 0;
          for (int m=0; m<k; ++m) {
            space[jj] += U[n_local * m + j] * Ainv[jj + k * m];
          }
        } // ends for over jj

        double d = 1.0;
        for (int jj=0; jj<k; ++jj) {
          d += space[jj] * U[n_local * jj + j];
        }
        for (int ll=0; ll<k; ++ll) { space[ll] = U[n_local * ll + j]; }
        char   TRANS = 'N';
        double ONE   = 1.0, ZERO = 0.0;
        int    INCX  = 1;
        dgemv_(&TRANS, &k, &k, &ONE, Ainv, &k, &U[j], &n_local, &ZERO, &space[0], &INCX);

        std::shared_ptr<double> space2 (new double[k]);
        for (int pp=0; pp<k; ++pp) { space2.get()[pp] = space[pp]; }
        dgemv_(&TRANS, &k, &k, &ONE, Ainv, &k, &(space2.get()[0]), &INCX, &ZERO,
                                               &space[0],  &INCX); 

        double adj = 0;
        for (int ii=0; ii<k; ++ii) { adj += U[j + n_local * ii] * space[ii]; }
        adj /= d;
 
        // test column j
        if (trhinv-adj <= trinv) {
          for (int prc=0; prc<size; ++prc) {
            // send data
            for (int i=0; i<k; ++i) u_vec[i] = U[j + n_local * i];
            u_vec[k] = j_global;
            MPI_Send(&u_vec[0], k + 1, MPI_DOUBLE, prc, 1, MPI_COMM_WORLD);
          }
        }
      } // ends not already selected if statement

      // test if a column has been found
      MPI_Test(&request, &done, &status);

    } // end of search loop

    // make sure cores are all working on the same iteration
    MPI_Barrier(MPI_COMM_WORLD);

    // send column (row in this formulation of U) to update A
    int find_root = fcol * size / n;

    // update A
    if (root) {

      // only record if this iteration actually found a new entry
      if (record) {
        for (int ii=0; ii<k; ++ii) {
          for (int y=0; y<k; ++y) {
            A[k * y + ii] += u_vec[ii] * u_vec[y];
          }
        }
        pi[int(u_vec[k])] = r + 1;
      }
    }
  } // end of iteration for loop

  // column selection is complete

  if (root) {
    // record selection
    if (VERBOSE) {
      std::cout << "pi: " << std::endl;
      for (int i=0; i<n; ++i) { std::cout << pi[i] << std::endl; }
      std::cout << std::endl;
    }

    std::shared_ptr<int> pi2 (new int[l]);
    convertPi(&pi[0], pi2, n, l);
    if (VERBOSE) {
      std::cout << "pi2: " << std::endl;
      for (int i=0; i<l; ++i) {std::cout << pi2.get()[i] << std::endl; }
      std::cout << std::endl;
    }
    if (WRITE) {
      std::ofstream pfile;
      pfile.open("p.txt");
      for (int i=0; i<l; ++i) {
        pfile << pi2.get()[i] << "\n";
      }
      pfile.close();
    }
  }

  std::cout << "Processor " << rank << " has finished" << std::endl;

  // end timer here
  if (root) {
    seconds = read_timer() - seconds;
    std::cout << "Computation time: " << seconds << " seconds on " << size << " core(s)" << std::endl;
  }

  MPI_Finalize();
  return 0;
}

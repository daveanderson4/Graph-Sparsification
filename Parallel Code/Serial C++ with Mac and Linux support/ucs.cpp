#include "sparsifier.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <float.h>
#include <stdlib.h>
#include <cassert>

using namespace std;

void convertPi(int* pi_in, int* pi_out, int n, int k) {
  for (int i=0; i<n; ++i) {
    if (pi_in[i] != 0) {
      pi_out[pi_in[i]-1] = i+1;
    }
  }
}

int main(int argc, char* argv[]) {

  // parameters
  int n = 73;
  int k = 35;
  int l = 2*k;
  if (argc > 3) {
    n = atoi(argv[1]);
    k = atoi(argv[2]);
    l = atoi(argv[3]);
  }
  bool VERBOSE = false;
  bool WRITE   = true;

  // read in data
  ifstream fin;
  fin.open ("U.txt");
  assert(fin);
  double *U = new double[n*k];
  for (int i=0; i<n; ++i) {
    for (int j=0; j<k; ++j) {
      fin >> U[i+n*j];
    }
  }
  fin.close();

  // run algorithm
  sparsifier S(n,k,l);
  int* pi = new int[n];
  for (int i=0; i<n; ++i) { pi[i] = 0; }
  S.columnSelect(pi, U);

  // record selection
  if (VERBOSE) {
    cout << "pi: " << endl;
    for (int i=0; i<n; ++i) { cout << pi[i] << endl; }
    cout << endl;
  }
  int* pi2 = new int[l];
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
  delete U;
  delete pi;
  delete pi2;
  return 0;
}

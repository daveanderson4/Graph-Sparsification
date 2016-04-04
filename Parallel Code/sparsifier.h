#ifndef _SPARSIFIER_H_
#define _SPARSIFIER_H_

#include <vector>
#include <memory>

class sparsifier {
  public:
    // constructor
    // problem dimensions are set and T is calculated
    sparsifier(int n, int k, int l);

    // destructor
    ~sparsifier();

    // perform column selection algorithm on U
    void columnSelect(int* pi, double* U);

    // return class member variable n
    int getN();

    // return class member variable T
    double getT();

    // lambda is the root of this function
    static double root1(double x, std::vector<double>& s, int n, int k, double T, int r, double c);

    // lambda_hat is the root of this function
    static double root2(double lh, std::vector<double>& s, int n, int k, double lambda, int r, double c);

    // a constant for root2
    static double root2helper(std::vector<double>& s, int n, int k, double lambda, int r);

    // for the inverse of a small (k-by-k) matrix
    static void inverse(std::shared_ptr<double> A, int N);

  private:
    // calculates problem parameter T
    double findT(int n, int k, int l);

    // member variables
    int n;
    int k;
    int l;
    double T;
};
#endif


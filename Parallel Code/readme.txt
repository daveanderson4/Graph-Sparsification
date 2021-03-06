A Parallelized Graph Sparsification Example
David G. Anderson
2016

Autonomous Systems Example
Data from http://snap.stanford.edu/data/index.html

In this version of the code, step 3 of the unweighted
graph sparsification algorithm is performed in parallel,
which is the bottleneck of computation.  Step 3 takes
the computed lambda_hat and searches for the next usable
column.

Computations were run on multiple cores of a single node
of NERSC’s Cori.  The problem has dimensions:
n=1145, k=492, l=493

Data files:
 - U.txt
 - as.txt

Computation files:
 - makefile
 - ucs.cpp
     main algorithm, parallelized with MPI
 - sparsifier.cpp/h
     some support routines
 - solver.cpp/h
     root finders

Output files:
 - p.txt
     the selected columns
 - results_tree_plot.png
     plot of the spanning tree generated by the algorithm
 - results_time.png
     time comparison of the same problem run on differing
     numbers of cores
 - results_strong_scaling.png
     strong scaling efficiency
 - results_weak_scaling.png
     weak scaling efficiency

Other files:
 - AS.m
     a Matlab file performing preprocessing of orthogonal
     decomposition of the graph incidence matrix
 - as_graphs.R
     plotting support in R
 - as19981229.txt
     original data file

Other folders:
 - Matlab code
     Matlab implementation
 - Serial C++ with Mac and Linux support
     An earlier version of the code in serial form, which
     includes compilation support for Max in addition to 
     Linux
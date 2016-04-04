Code for an example with Autonomous Systems
Data from http://snap.stanford.edu/data/index.html

David G. Anderson
2016

files:
	0. Data
		as19981229.txt  - files containing graph edges
		as.txt
	1. Matlab
		AS.m 		- load raw data and preprocessing
		U.txt  		- output containing orthonormal matrix (see paper)
	2. C++
		makefile	- compilation for Cray and Apple machines provided
				  change compiler and flag options to recompile
		Sparsifier	- executable created from compilation
				- reads in U.txt from matlab code, outputs p.txt
				  with selected edges
				- required inputs: n, k, l
					- n is number of edges in starting graph
					- k is number of nodes - 1
					- l is desired number of edges (k < l <= n)
					- use ‘./Sparsifier 1145 492 493’ for AS example
		ucs.cpp 	- main function
		sparsifier.cpp	- column selection algorithm and subroutines
		sparsifier.h
		solver.cpp	- root finders
		solver.h
		p.txt		- output containing selected edges
	3. R
		as_graphs.R	- plots sparsifier over original graph using a 
				  force-directed graph plotting algorithm




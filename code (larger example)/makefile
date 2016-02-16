CSRC := ucs.cpp sparsifier.cpp solver.cpp
FSRC := 

CCRAY := CC
FCRAY := ftn
CLOCAL := g++
FLOCAL := ftn

CFLAGSCRAY := -pg -g -dynamic -O3 -openmp -std=c++11 -mkl=sequential
FFLAGSCRAY := -pg -g -dynamic -O3 -openmp -mkl=sequential
CFLAGSLOCAL := -Wall -llapack -lblas -O3 -std=c++11
FFLAGSLOCAL :=

CXX := $(CLOCAL)
FC := $(FLOCAL)
CFLAGS := $(CFLAGSLOCAL)
FFLAGS := $(FFLAGSLOCAL)

COBJ := $(CSRC:%.cpp=%.o)
FOBJ := $(FSRC:%.f=%.o)

TARGETS := Sparsifier

all: $(TARGETS)

$(COBJ) : %.o : %.cpp
	$(CXX) -c $< -o $@ $(CFLAGS)

$(FOBJ) : %.o : %.f
	$(FC) -c $< -o $@ $(FFLAGS)

Sparsifier: $(COBJ) $(FOBJ)
	$(CXX) $^ -o $@ $(CFLAGS)

clean:
	rm *.o



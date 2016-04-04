CSRC := ucs.cpp sparsifier.cpp solver.cpp
FSRC := 

CXX := CC
FC := ftn

CFLAGS := -pg -g -dynamic -O3 -openmp -std=c++11 -mkl
FFLAGS := -pg -g -dynamic -O3 -openmp -mkl

COBJ := $(CSRC:%.cpp=%.o)
FOBJ := $(FSRC:%.f=%.o)
CTESTOBJ := $(CTEST:%.cpp=%.o)

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



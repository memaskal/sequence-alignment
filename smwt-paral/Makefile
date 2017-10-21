CCX=g++
CXXFLAGS=-Wall -Wextra -O3
OMPFLAGS=-fopenmp
EXTRAFLAGS=-march=native -mfpmath=sse

all: smwt-paral.cpp
	$(CXX) -o smwt-paral $(CXXFLAGS) $(OMPFLAGS) $(EXTRAFLAGS) smwt-paral.cpp



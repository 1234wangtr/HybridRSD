# Hybrid approach for Regular Syndrome Decoding Problem

The C++ project is a simple demo with (f,u)-guess(guessing fu error-free positions in f blocks) for our hybrid approach for RSD.
For the sake of simplicity, we use the normal version of Gaussian elimination with a complexity of $O(n^3)$ instead of $O(n^3/logn)$.
The Python script for estimating complexity is in "rsd\_newmodel.py" and the quick version is in "rsd\_quick\_esti.py", where the apis are hybrid\_bigq\_quick and hybrid\_2\_quick.

## Parameter

* RSD\_N,RSD\_K,RSD\_H,RSD\_B are parameters for RSD. F,U,Q are parameters f,u,g for our algorithm.

## Run

* g++ -O3 -o test.exe par\_matrix.cpp m\_matrix.cpp full\_algeb\_attack.cpp

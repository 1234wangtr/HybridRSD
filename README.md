# Hybrid approach for Regular Syndrome Decoding Problem

The C++ project is a simple demo with (f,u)-guess(guessing fu error-free positions in f blocks) for our hybrid approach for binary RSD. 
For the sake of simplicity, we use the normal version of Gaussian elimination with a complexity of $O(n^3)$ instead of $O(n^3/logn)$.
The Python script for estimating complexity is in "rsd_newmodel.py".

## Parameter

- RSD_N,RSD_K,RSD_H,RSD_B are parameters for RSD. F,U,Q are parameters f,u,g for our algorithm.

## Run

- g++ -O3 -o test.exe par_matrix.cpp m_matrix.cpp full_algeb_attack.cpp


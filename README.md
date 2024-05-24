# Hybrid approach for Regular Syndrome Decoding Problem

This is a simple demo without further optimizations for our hybrid approach for binary RSD. 
For the sake of simplicity, we use the normal version of Gaussian elimination with a complexity of $O(n^3)$ instead of $O(n^3/logn)$.

## Parameter

- RSD_N,RSD_K,RSD_H,RSD_B are parameters for RSD. F,U,Q are parameters for our algorithm.

## Run

- g++ -O3 -o test.exe par_matrix.cpp m_matrix.cpp full_algeb_attack.cpp


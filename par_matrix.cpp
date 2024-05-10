#include "par_matrix.h"

#include <iostream>
using namespace std;

int dimTrans(int bg, int ed, int n) {
	if (bg > ed) {
		swap(bg, ed);
	}
	int sum = bg * (n  + n - bg+1) / 2;
	sum += ed-bg;
	return sum;
}

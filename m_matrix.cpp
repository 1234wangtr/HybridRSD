#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <vector>
#include <algorithm> // std::move_backward
#include <random> // std::default_random_engine
#include <chrono> // std::chrono::system_clock

#include <map>
#include <assert.h>
using namespace std;

#include "m_matrix.h"

//const int n = 1280;
//const int k = 640;
//const int r = 640;

const int n = 100;
const int k = 50;
const int r = 50;

time_t bg, ed;


// 尝试浪费空间的
//char matrix[k][n];
char** matrix;

Matrix HWhole(k, n);
Matrix HCopy(k,n);
Matrix GWhole(r, n);
Matrix GCopy(r,n);

void init() {
	matrix = new char*[k];
	for (int i = 0; i < k; i++) {
		matrix[i] = new char[n];
	}
}

void clearAll() {
	for (int i = 0; i < k; i++) {
		delete [] matrix[i];
	}
	delete[] matrix;
}

Vec vecAdd(const Vec& v1, const Vec& v2) {
	assert(v1.size == v2.size);
	Vec ret(v1.size);
	for (int i = 0; i < v1.size; i++) {
		ret.data[i] = v1.data[i] + v2.data[i];
	}
	return ret;
}

Vec matrixMult(const Matrix & m, const Vec& v) {

	int h = m.height;
	int w = m.width;
	//printf("h=%d w=%d\n", h, w);
	Vec ret(h);
	for (int i = 0; i < h; i++) {
		int res = 0;
		for (int j = 0; j < w; j++) {
			res += m.data[i][j] * v.data[j];
		}
		ret.data[i] = res % 2;
	}
	return ret;
}

// 需要传入一个向量 以记录打乱次序
Matrix randomGauss(Matrix mo, int l, vector<int>& v) {
	Matrix m(mo.height, mo.width);
	//random
	//std::vector<int> v;
	v.clear();
	int width = mo.width;
	int height = mo.height;
	for (int i = 0; i < width; ++i) {
		v.push_back(i);
	}
	// obtain a time-based seed:
	
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	seed = rand(); // get rid of time 
	std::shuffle(v.begin(), v.end(), std::default_random_engine(seed));
	

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			m.data[j][i] = mo.data[j][v[i]];
		}
	}
	//debug
	//printf("after random:\n");
	//m.printIt();


	//左上角l-r gauss消元
	int leftUp = height - l;
	for (int i = 0; i < leftUp; i++) {
		//printf("i=%d\n", i);
		//m.printIt();

		//find 1
		//cout << "dbg mdata=" << (int)m.data[i][i] << endl;
		if (m.data[i][i] == 0) {
			int find = -1;
			for (int j = i + 1; j < height; j++) {
				if (m.data[j][i] == 1) {
					find = j;
					break;
				}
			}
			if (find == -1) {
				//GG
				//cout << "GGGG" << endl;
				return Matrix(1, 1);
			}


			//add to line i
			for (int j = i; j < width; j++) {
				m.data[i][j] = (m.data[i][j] + m.data[find][j]) % 2;
			}
			//cout << "swap over" << endl;
		}
		//printf("after add\n");
		//m.printIt();
		for (int j = 0; j < height; j++) {
			if (j == i)continue;
			if (m.data[j][i] == 1) {
				for (int k = 0; k < width; k++) {
					m.data[j][k] = (m.data[j][k] + m.data[i][k]) % 2;
				}
			}
		}
	}

	//debug
	/*
	printf("after gaussaaaaaaaaaaaaaa:\n");
	m.printIt();
	printf(".......");*/
	return m;
}

Vec getVec(int len, int weight) {
	//generate a vec with len 0/1 and weight 1
	vector<int> v;
	v.clear();
	for (int i = 0; i < len; ++i) {
		v.push_back(i);
	}
	// obtain a time-based seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle(v.begin(), v.end(), std::default_random_engine(seed));

	Vec vv(len);
	for (int i = 0; i < len; i++) {
		vv.data[i] = 0;
	}
	for (int i = 0; i < weight; i++) {
		vv.data[v[i]] = 1;
	}
	return vv;
}

void printMarix() {
	ofstream output("dbg.txt", ios::trunc);

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < n; j++) {
			output << (int)(matrix[i][j]) << " ";
		}
		output << endl;
	}
}

Matrix getSub(Matrix& m, int bgh, int bgw, int hei, int wid) {
	Matrix mm(hei, wid);
	if (hei + bgh > m.height || wid + bgw > m.width) {
		printf("ERR:bgh=%d hei=%d bgw=%d wid=%d m.hei=%d m.wei=%d\n", bgh, hei, bgw, wid, m.height, m.width);
		assert(false);
	}
	for (int i = 0; i < hei; i++) {
		for (int j = 0; j < wid; j++) {
			mm.data[i][j] = m.data[i + bgh][j + bgw];
		}
	}
	return mm;
}


Matrix getTrans(Matrix & m) {
	Matrix mm(m.width, m.height);
	for (int i = 0; i < m.width; i++) {
		for (int j = 0; j < m.height; j++) {
			mm.data[i][j] = m.data[j][i];
		}
	}
	return mm;
}


Matrix matMulMat(const Matrix & m, const Matrix & mm) {
	assert(m.width == mm.height);
	Matrix res(m.height, mm.width);
	for (int i = 0; i < m.height; i++) {
		for (int j = 0; j < mm.width; j++) {
			res.data[i][j] = 0;
		}
	}
	for (int i = 0; i < m.height; i++) {
		for (int j = 0; j < m.width; j++) {
			for (int k = 0; k < mm.width; k++) {
				res.data[i][k] ^= m.data[i][j] * mm.data[j][k];
			}
		}
	}
	return res;
}

Matrix justGauss(Matrix mo) {
	Matrix m(mo.height, mo.width);
	int width = mo.width;
	int height = mo.height;
	
	
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			m.data[j][i] = mo.data[j][i];
		}
	}
	//debug
	//printf("after random:\n");
	//m.printIt();


	//左上角l-r gauss消元
	int leftUp = height ;

	// another type of matrix
	//if (height > width)leftUp = width;


	bool flag = true;
	int fail_at = -1;

	for (int i = 0; i < leftUp; i++) {
		//printf("i=%d\n", i);
		//m.printIt();

		//find 1
		//cout << "dbg mdata=" << (int)m.data[i][i] << endl;
		if (m.data[i][i] == 0) {
			int find = -1;
			for (int j = i + 1; j < height; j++) {
				if (m.data[j][i] == 1) {
					find = j;
					break;
				}
			}
			if (find == -1) {
				//GG
				//cout << "GGGG when i="<<i <<" goal="<<leftUp<< endl;
				return Matrix(1, 1);
				flag = false;
				fail_at = i;
				justGaussDbg(mo, fail_at);
				return Matrix(1, 1);

				//
				for (int iii = 0; iii < i + 8; iii++) {
					for (int jjj = 0; jjj < i+8; jjj++) {
						cout << (int)(mo.data[iii][jjj]);
					}cout << endl;
				}cout << endl;
				cout << "------------" << endl;
				for (int iii = 0; iii < i+8; iii++) {
					for (int jjj = 0; jjj < i+8; jjj++) {
						cout << (int)(m.data[iii][jjj]);
					}cout << endl;
				}cout << endl;
				return Matrix(1, 1);
			}


			//add to line i
			for (int j = i; j < width; j++) {
				m.data[i][j] = (m.data[i][j] + m.data[find][j]) % 2;
			}
			//cout << "swap over" << endl;
		}
		//printf("after add\n");
		//m.printIt();
		for (int j = 0; j < height; j++) {
			if (j == i)continue;
			if (m.data[j][i] == 1) {
				for (int k = 0; k < width; k++) {
					m.data[j][k] = (m.data[j][k] + m.data[i][k]) % 2;
				}
			}
		}
	}

	return m;
}

Matrix justGaussDbg(Matrix mo,int idx) {
	Matrix m(mo.height, mo.width);
	int width = mo.width;
	int height = mo.height;


	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			m.data[j][i] = mo.data[j][i];
		}
	}
	//debug
	//printf("after random:\n");
	//m.printIt();


	//左上角l-r gauss消元
	int leftUp = height;

	// another type of matrix
	//if (height > width)leftUp = width;


	bool flag = true;
	int fail_at = -1;

	for (int i = 0; i < leftUp; i++) {
		//printf("i=%d\n", i);
		//m.printIt();

		//find 1
		//cout << "dbg mdata=" << (int)m.data[i][i] << endl;
		//dbg
		if (i + 5 >= idx) {
			for (int iii = 0; iii < i + 8; iii++) {
				for (int jjj = 0; jjj < i + 8; jjj++) {
					cout << (int)(mo.data[iii][jjj]);
				}cout << endl;
			}cout << endl;
			cout << "------------" << endl;
			for (int iii = 0; iii < i + 8; iii++) {
				for (int jjj = 0; jjj < i + 8; jjj++) {
					cout << (int)(m.data[iii][jjj]);
				}cout << endl;
			}cout << endl;
		}
		if (m.data[i][i] == 0) {
			int find = -1;
			for (int j = i + 1; j < height; j++) {
				if (m.data[j][i] == 1) {
					find = j;
					break;
				}
			}
			if (find == -1) {
				//GG
				//cout << "GGGG when i=" << i << " goal=" << leftUp << endl;
				
				return Matrix(1, 1);
			}


			//add to line i
			for (int j = i; j < width; j++) {
				m.data[i][j] = (m.data[i][j] + m.data[find][j]) % 2;
			}
			//cout << "swap over" << endl;
		}
		//printf("after add\n");
		//m.printIt();
		for (int j = 0; j < height; j++) {
			if (j == i)continue;
			if (m.data[j][i] == 1) {
				for (int k = 0; k < width; k++) {
					m.data[j][k] = (m.data[j][k] + m.data[i][k]) % 2;
				}
			}
		}
	}

	return m;
}



Matrix figureEqus(Matrix mo, Vec& v) {
	Matrix m(mo.height, mo.width);
	int width = mo.width;
	int height = mo.height;


	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			m.data[j][i] = mo.data[j][i];
		}
	}
	//debug
	//printf("after random:\n");
	//m.printIt();


	//左上角l-r gauss消元
	int leftUp = max(height,width);
	for (int i = 0; i < leftUp; i++) {
		//printf("i=%d\n", i);
		//m.printIt();

		//find 1
		//cout << "dbg mdata=" << (int)m.data[i][i] << endl;
		if (m.data[i][i] == 0) {
			int find = -1;
			for (int j = i + 1; j < height; j++) {
				if (m.data[j][i] == 1) {
					find = j;
					break;
				}
			}
			if (find == -1) {
				//GG
				//cout << "GGGG" << endl;
				return Matrix(1, 1);
			}


			//add to line i
			for (int j = i; j < width; j++) {
				m.data[i][j] = (m.data[i][j] + m.data[find][j]) % 2;
			}
			v.data[i] ^= v.data[find];
			//cout << "swap over" << endl;
		}
		//printf("after add\n");
		//m.printIt();
		for (int j = 0; j < height; j++) {
			if (j == i)continue;
			if (m.data[j][i] == 1) {
				for (int k = 0; k < width; k++) {
					m.data[j][k] = (m.data[j][k] + m.data[i][k]) % 2;
				}
				v.data[j] ^= v.data[i];
			}
		}
	}

	return m;
}


int vecDot(const Vec& v1, const Vec& v2) {
	assert(v1.size == v2.size);
	int ret = 0;
	for (int i = 0; i < v1.size; i++) {
		ret += v1.data[i] * v2.data[i];
	}
	return ret % 2;
}



Matrix randomGauss_my(Matrix mo, int l, vector<int>& v) {
	//printf("random Gauss\n");
	Matrix m(mo.height, mo.width);
	//random
	v.clear();
	int width = mo.width;
	int height = mo.height;
	for (int i = 0; i < width; ++i) {
		v.push_back(i);
	}
	// obtain a time-based seed:

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	

	// 作弊！！！
	/*swap(v[0], v[n - 1]);
	swap(v[1], v[n - 2]);
	swap(v[2], v[n - 3]);
	swap(v[3], v[n - 4]);
	swap(v[4], v[n - 5]);
	swap(v[5], v[n - 6]);
	swap(v[6], v[n - 7]);
	swap(v[7], v[n - 8]);
	swap(v[8], v[n - 9]);
	swap(v[9], v[n - 10]);*/

	int chonghe = 11;
	auto bg = v.begin();
	auto ed = v.end();
	/*
	for (int i = 0; i < chonghe; i++) {
		bg++;
		ed--;
		swap(v[i], v[n - i - 1]);
	}
	bg++; ed--;
	//bg++; bg++;
	//ed--; ed--;

	//swap(v[n-11-1], v[92]);
	//swap(v[11], v[23]);    // very important !
	*/
	seed = rand(); // get rid of time 
	std::shuffle(bg,ed, std::default_random_engine(seed));
	//printf("v %d=%d   v %d=%d\n", 11, v[11], n - 11 - 1, v[n - 11 - 1]);
	/*
	for (int i = 0; i < n; i++) {
		printf("%d ", v[i]);
	}cout << endl;*/

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			m.data[j][i] = mo.data[j][v[i]];
		}
	}
	//debug
	//printf("after random:\n");
	//m.printIt();


	//左上角l-r gauss消元
	int leftUp = height - l;
	for (int i = 0; i < leftUp; i++) {
		//printf("i=%d\n", i);
		//m.printIt();

		//find 1
		//cout << "dbg mdata=" << (int)m.data[i][i] << endl;
		if (m.data[i][i] == 0) {
			int find = -1;
			for (int j = i + 1; j < height; j++) {
				if (m.data[j][i] == 1) {
					find = j;
					break;
				}
			}
			if (find == -1) {
				//GG
				//cout << "GGGG" << endl;
				return Matrix(1, 1);
			}


			//add to line i
			for (int j = i; j < width; j++) {
				m.data[i][j] = (m.data[i][j] + m.data[find][j]) % 2;
			}
			//cout << "swap over" << endl;
		}
		//printf("after add\n");
		//m.printIt();
		for (int j = 0; j < height; j++) {
			if (j == i)continue;
			if (m.data[j][i] == 1) {
				for (int k = 0; k < width; k++) {
					m.data[j][k] = (m.data[j][k] + m.data[i][k]) % 2;
				}
			}
		}
	}

	//debug
	/*
	printf("after gaussaaaaaaaaaaaaaa:\n");
	m.printIt();
	printf(".......");*/
	return m;
}


Matrix justGauss2(Matrix mo) {
	Matrix m(mo.height, mo.width);
	int width = mo.width;
	int height = mo.height;


	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			m.data[j][i] = mo.data[j][i];
		}
	}
	//debug
	//printf("after random:\n");
	//m.printIt();


	//左上角l-r gauss消元
	int leftUp = height;

	// another type of matrix
	if(leftUp > width - 1)
		leftUp = width-1;

	for (int i = 0; i < leftUp; i++) {
		//printf("i=%d\n", i);
		//m.printIt();

		//find 1
		//cout << "dbg mdata=" << (int)m.data[i][i] << endl;
		if (m.data[i][i] == 0) {
			int find = -1;
			for (int j = i + 1; j < height; j++) {
				if (m.data[j][i] == 1) {
					find = j;
					break;
				}
			}
			if (find == -1) {
				//GG
				//cout << "GGGG at "<<i <<" goal="<<leftUp<< endl;
				return Matrix(1, 1);
			}


			//add to line i
			for (int j = i; j < width; j++) {
				m.data[i][j] = (m.data[i][j] + m.data[find][j]) % 2;
			}
			//cout << "swap over" << endl;
		}
		//printf("after add\n");
		//m.printIt();
		for (int j = 0; j < height; j++) {
			if (j == i)continue;
			if (m.data[j][i] == 1) {
				for (int k = 0; k < width; k++) {
					m.data[j][k] = (m.data[j][k] + m.data[i][k]) % 2;
				}
			}
		}
	}

	return m;
}

Matrix justGaussLeftUp(Matrix mo,int lu) {
	Matrix m(mo.height, mo.width);
	int width = mo.width;
	int height = mo.height;


	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			m.data[j][i] = mo.data[j][i];
		}
	}
	//debug
	//printf("after random:\n");
	//m.printIt();


	//左上角l-r gauss消元
	int leftUp = lu;

	// another type of matrix
	if (leftUp > width - 1)
		leftUp = width - 1;

	for (int i = 0; i < leftUp; i++) {
		//printf("i=%d\n", i);
		//m.printIt();

		//find 1
		//cout << "dbg mdata=" << (int)m.data[i][i] << endl;
		if (m.data[i][i] == 0) {
			int find = -1;
			for (int j = i + 1; j < height; j++) {
				if (m.data[j][i] == 1) {
					find = j;
					break;
				}
			}
			if (find == -1) {
				//GG
				//cout << "GGGG at " << i << " goal=" << leftUp << endl;
				return Matrix(1, 1);
			}


			//add to line i
			for (int j = i; j < width; j++) {
				m.data[i][j] = (m.data[i][j] + m.data[find][j]) % 2;
			}
			//cout << "swap over" << endl;
		}
		//printf("after add\n");
		//m.printIt();
		for (int j = 0; j < height; j++) {
			if (j == i)continue;
			if (m.data[j][i] == 1) {
				for (int k = 0; k < width; k++) {
					m.data[j][k] = (m.data[j][k] + m.data[i][k]) % 2;
				}
			}
		}
	}

	return m;
}



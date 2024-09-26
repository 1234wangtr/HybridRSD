#pragma once
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

/*
	The m_matrix.h as well as the m_matrix.cpp are about matrix, vector and Gaussian elimination
*/


class Matrix {
public:
	int width;
	int height;
	char** data;
	Matrix(int _height = 10, int _width = 10) {
		height = _height;
		width = _width;

		data = new char*[height];
		for (int i = 0; i < height; i++) {
			data[i] = new char[width];
		}
	}
	~Matrix() {
		for (int i = 0; i < height; i++) {
			delete[] data[i];
		}
		delete[]data;

	}
	Matrix(const Matrix & m) {
		height = m.height;
		width = m.width;


		data = new char*[height];
		for (int i = 0; i < height; i++) {
			data[i] = new char[width];
			for (int j = 0; j < width; j++) {
				data[i][j] = m.data[i][j];
			}
		}
	}

	Matrix& operator=(const Matrix& m) {
		for (int i = 0; i < height; i++) {
			delete[] data[i];
		}
		delete[]data;

		height = m.height;
		width = m.width;


		data = new char*[height];
		for (int i = 0; i < height; i++) {
			data[i] = new char[width];
			for (int j = 0; j < width; j++) {
				data[i][j] = m.data[i][j];
			}
		}
		return *this;
	}

	bool isEmpty() {
		if (width*height > 1)return false;
		return true;
	}
	void randomIt() {
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				data[i][j] = rand() % 2;
			}
		}
	}

	void printIt() {
		printf("width=%d height=%d\n", width, height);
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				printf("%d ", data[i][j]);
			}
			printf("\n");
		}
	}



};

Matrix getSub(Matrix& m, int bgh, int bgw, int hei, int wid);

Matrix getTrans(Matrix & m);

class Vec {
public:
	int size;
	int* data;
	Vec(string s) {
		size = s.length();
		data = new int[size];
		for (int i = 0; i < size; i++) {
			data[i] = s[i] - '0';
		}
	}

	Vec(int s) {
		size = s;
		data = new int[size];
	}
	~Vec() {
		delete[] data;
	}
	Vec(const Vec& v) {
		size = v.size;
		data = new int[size];
		for (int i = 0; i < size; i++) {
			data[i] = v.data[i];
		}
	}

	Vec& operator=(const Vec& v) {
		delete[] data;
		size = v.size;
		data = new int[size];
		for (int i = 0; i < size; i++) {
			data[i] = v.data[i];
		}
		return *this;
	}

	Vec operator +(const Vec& v) {
		assert(size == v.size);
		for (int i = 0; i < size; i++) {
			data[i] += v.data[i];
		}
		return *this;
	}

	void printIt() const {
		cout << "vec size=" << size << endl;
		for (int i = 0; i < size; i++) {
			cout << (int)data[i];
		}cout << endl;
	}

	void printIt(const int cycle) const {
		cout << "vec size=" << size << endl;
		for (int i = 0; i < size; i++) {
			cout << (int)data[i];
			if (i % cycle == cycle - 1)cout << " ";
		}cout << endl;
	}

	void printItOri() {
		cout << "vec size=" << size << endl;
		for (int i = 0; i < size; i++) {
			cout << (int)data[i] << " ";
		}
	}

	bool operator <(const Vec& ano) const {
		for (int i = 0; i < size; i++) {
			if (data[i] < ano.data[i]) {
				return true;
			}
			if (data[i] > ano.data[i]) {
				return false;
			}
		}
		return false;
	}

	bool operator == (const Vec& ano) const {
		if (size == ano.size) {
			for (int i = 0; i < size; i++) {
				if (data[i] != ano.data[i])return false;
			}
			return true;
		}
		else {
			return false;
		}
	}

	int getWeight() {
		int res = 0;
		for (int i = 0; i < size; i++) {
			if (data[i] == 1)res++;
		}
		return res;
	}
};

void init();

void clearAll();


Vec vecAdd(const Vec& v1, const Vec& v2);

int vecDot(const Vec& v1, const Vec& v2);

Matrix randomGauss(Matrix mo, int l, vector<int>& v);
Matrix randomGauss_my(Matrix mo, int l, vector<int>& v);

Matrix figureEqus(Matrix,Vec&); 
Matrix justGauss(Matrix);
Matrix justGauss2(Matrix);
Matrix justGaussDbg(Matrix,int);

Matrix justGaussLeftUp(Matrix mo, int lu);

Vec getVec(int len, int weight);

void printMarix();
Vec matrixMult(const Matrix & m, const Vec& v);

Matrix matMulMat(const Matrix & m, const Matrix & mm);

extern string tmpAns;
extern time_t bg, ed;

extern const int n;
extern const int k;
extern const int r;

extern char** matrix;

extern Matrix HWhole;
extern Matrix HCopy;
extern Matrix GWhole;
extern Matrix GCopy;

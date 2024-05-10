#pragma once

#include <iostream>
#include "m_matrix.h"
using namespace std;
int dimTrans(int bg, int ed, int tot);

class ParMat {
public:
	int row=100;
	int varNum;
	int parNum;	//1 is also a par
	char** var_var;
	char** par_var;
	char** par_par;

	Matrix* gaussedParM;



	ParMat(int _row, int _varNum, int _parNum) {
		row = _row;
		varNum = _varNum;
		parNum = _parNum;
		var_var = new char* [row];	//v0v0 v0v1 ... v0vn ;  v1v2 ... v1vn ... 
		par_var = new char* [row];
		par_par = new char* [row];

		for (int i = 0; i < row; i++) {
			var_var[i] = new char[(varNum*varNum + varNum)/2];
			for (int j = 0; j < (varNum * varNum + varNum) / 2; j++)var_var[i][j] = 0;
			
			par_var[i] = new char[parNum * varNum];
			for (int j = 0; j < parNum * varNum; j++)par_var[i][j] = 0;

			par_par[i] = new char[(parNum * parNum + parNum) / 2];
			for (int j = 0; j < (parNum * parNum + parNum) / 2; j++)par_par[i][j] = 0;
		}
		gaussedParM = nullptr;
	}

	~ParMat() {
		for (int i = 0; i < row; i++) {
			delete [] var_var[i];	
			delete [] par_var[i];
			delete [] par_par[i];
		}
		delete[] var_var;
		delete[] par_var;
		delete[] par_par;
		if (gaussedParM) delete gaussedParM;
	}

	void addRow(int idx, char* var1, char* par1, char* var2, char* par2) {
		//init the idx row

		for (int i = 0; i < varNum; i++) {
			for (int j = 0; j < varNum; j++) {
				int dim = dimTrans(i, j, varNum);
				var_var[idx][dim] ^= var1[i] * var2[j];
			}
		}

		for (int i = 0; i < varNum; i++) {
			for (int j = 0; j < parNum; j++) {
				par_var[idx][i * parNum + j] ^= var1[i] * par2[j];
				par_var[idx][i * parNum + j] ^= var2[i] * par1[j];
			}
		}

		for (int i = 0; i < parNum; i++) {
			for (int j = 0; j < parNum; j++) {
				int dim = dimTrans(i, j, parNum);
				par_par[idx][dim] ^= par1[i] * par2[j];
				//par_par[idx][dim] ^= par1[j] * par2[i];
			}
		}

		//fresh    in par x var, 1*vi = vi*vi, so put it in var x var
		// par[-1] correspond to constant 1
		//for (int i = 0; i < varNum; i++) {
		//	var_var[idx][dimTrans(i, i, varNum)] ^= par_var[idx][i * parNum + parNum - 1];
		//}


		return;
	}

	void printRow(int idx) {
		cout << endl<< "row dbg:" << idx << endl;
		cout << "varvar:" << endl;
		for (int i = 0; i < varNum; i++) {
			for (int j = 0; j < varNum; j++) {
				int dim = dimTrans(i, j, varNum);
				cout << (int)var_var[idx][dim] << " ";
			}cout << endl;
		}

		cout << "parvar:" << endl;
		for (int i = 0; i < varNum; i++) {
			for (int j = 0; j < parNum; j++) {
				cout << (int)par_var[idx][i * parNum + j] << " ";
			}cout << endl;
		}

		cout << "parpar:" << endl;
		for (int i = 0; i < parNum; i++) {
			for (int j = 0; j < parNum; j++) {
				int dim = dimTrans(i, j, parNum);
				cout << (int)par_par[idx][dim] << " ";
			}cout << endl;
		}

	}

	int gaussPrepare() {
		//Matrix m(row, );
		int width = (varNum * varNum - varNum) / 2 + (parNum)*varNum + (parNum * parNum - parNum) / 2 + 1;
		int height = row;
		//TODO
		Matrix m(height, width);
		for (int i = 0; i < height; i++) {
			// var var(no v^2) || par var || par par(no p^2, p[-1]->1)
			int idx = 0;
			for (int j = 0; j < varNum; j++) {
				for (int l = j + 1; l < varNum; l++) {
					int dim = dimTrans(j, l, varNum);
					m.data[i][idx++] = var_var[i][dim];
				}
			}

			for (int j = 0; j < varNum; j++) {
				for (int l = 0; l < parNum - 1; l++) {
					m.data[i][idx++] = par_var[i][j * parNum + l];
				}
				m.data[i][idx++] = par_var[i][j * parNum + parNum - 1] ^ var_var[i][dimTrans(j,j,varNum)];
			}

			// add p^2 to 1*p (1*1 in the end)
			for (int j = 0; j < parNum; j++) {
				for (int l = j + 1; l < parNum; l++) {
					m.data[i][idx] = par_par[i][dimTrans(j, l, parNum)];
					if (l == parNum - 1)m.data[i][idx] ^= par_par[i][dimTrans(j, j, parNum)];
					idx++;
				}
			}
			m.data[i][idx++] = par_par[i][dimTrans(parNum - 1, parNum - 1, parNum)];
			assert(idx == width);
		}

		//cout << "dbg before gauss matrix:" << endl;
		//m.printIt();

		Matrix gaussedM = justGaussLeftUp(m,(varNum*varNum-varNum)/2 );
		//cout << "dbg after gauss matrix:" << endl;
		//gaussedM.printIt();

		int newH = height - (varNum * varNum - varNum) / 2;
		int newW = width - (varNum * varNum - varNum) / 2;
		gaussedParM = new Matrix(newH,newW);
		for (int i = 0; i < newH; i++) {
			for (int j = 0; j < newW; j++) {
				gaussedParM->data[i][j] = gaussedM.data[i + (varNum * varNum - varNum) / 2][j + (varNum * varNum - varNum) / 2];
			}
		}
		//gaussedParM->printIt();

		return 0;
	}



};
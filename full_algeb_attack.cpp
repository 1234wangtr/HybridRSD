#include "m_matrix.h"
#include "par_matrix.h"

using namespace std;


const int RSD_N = 400;
const int RSD_K = 139;
const int RSD_H = 25;
const int RSD_B = 16;

const int F = 20;
const int U = 1;

const int Q = 2;
const int U0 = Q*(RSD_B - 1);

Vec err(RSD_N + 1);
Vec err2(RSD_N - RSD_H + 1 - F * U);

Vec partAns(RSD_N);

int dbg = 0;
int important_print = 0;
bool rand_row = false;
bool permu_it = true;
int vrfy_dbg = 0;

bool fu_cheat = false;

int quadMult(const Vec & vbg, const Vec & ved, Vec& res) {
	int new_n = vbg.size ;
	int t = vbg.size - 1;
	int tmpIdx = 0;
	for (int i = 0; i < t - 1; i++) {
		for (int j = i + 1; j < t; j++) {
			int tmp = 0;
			tmp += vbg.data[i] * ved.data[j] + vbg.data[j] * ved.data[i];
			res.data[tmpIdx] = tmp % 2;
			tmpIdx++;
		}
	}
	for (int i = 0; i < t; i++) {
		int tmp = 0;
		tmp += vbg.data[i] * ved.data[i] + vbg.data[t] * ved.data[i] + vbg.data[i] * ved.data[t];
		res.data[tmpIdx] = tmp%2;
		tmpIdx++;
	}
	res.data[tmpIdx] = vbg.data[t] * ved.data[t];
	return 0;


	for (int jbg = 0; jbg < new_n; jbg++) {
		// degree = 1
		int idx_bg = jbg*(2 * new_n + 1 - jbg) / 2;
		cout << "jbg=" << jbg << " idxbg=" << idx_bg << endl;
		// e.g. 1 * a_t + a_t * 1 + a_t * a_t
		res.data[idx_bg] = vbg.data[new_n] * ved.data[jbg] + vbg.data[jbg] * ved.data[new_n] + vbg.data[jbg] * ved.data[jbg];
		res.data[idx_bg] %= 2;
		for (int jed = jbg + 1; jed < new_n; jed++) {
			idx_bg++;
			res.data[idx_bg] = vbg.data[jbg] * ved.data[jed] + vbg.data[jed] * ved.data[jbg];
			res.data[idx_bg] %= 2;
		}
	}
	return 0;
	
}

int guessParPos(int* pos, int tmpPos, Matrix* m) {
	if (tmpPos == Q) {
		//debug
		if (dbg) {
			cout << "dbg pos:" << endl;
			for (int i = 0; i < Q; i++) {
				cout << pos[i] << " ";
			}cout << endl;
		}

		// cheat
		bool cheat_ok = true;	// print if the guess is correct
		for (int i = 0; i < Q; i++) {
			if (err.data[(i + RSD_H - Q) * RSD_B + pos[i]] != 1) {
				cheat_ok = false;
				break;
			}
		}
		if (cheat_ok) {
			cout << "(cheat) This is the corrrect variance guess!!!" << endl;
		}
		else {
			//return -1;
		}
		
		// do sth
		Matrix mm(m->height, RSD_K - RSD_H - F * U - U0 + 1);
		for (int i = 0; i < m->height; i++) {
			// coeff of var
			for (int j = 0; j < RSD_K - RSD_H - F * U - U0; j++) {
				int sum = 0;
				//if(pos[j]<RSD_B-1)sum = m->data[i][j * (U0 + 1) + pos[j]];
				for (int l = 0; l < Q; l++) {
					if (pos[l] < RSD_B - 1) {
						sum += m->data[i][j * (U0 + 1) + l * (RSD_B - 1) + pos[l]];
					}
				}
				sum += m->data[i][j*(U0+1)+U0];
				sum %= 2;
				mm.data[i][j] = sum;
			}
			// constant
			int offset = (U0 + 1) * (RSD_K - RSD_H - F * U - U0);
			//cout << "offset=" << offset << endl;
			//cout << "m width=" << m->width << endl;
			
			int sum = 0;
			int vec[U0+1];
			for (int j = 0; j < U0 + 1; j++)vec[j] = 0;
			vec[U0] = 1;
			for (int j = 0; j < Q; j++)
				if(pos[j]<RSD_B-1)
					vec[(RSD_B - 1) * j + pos[j]] = 1;
			
			for (int j = 0; j < U0 + 1; j++) {
				for (int l = j+1; l < U0 + 1; l++) {
					sum += vec[j] * vec[l] * m->data[i][ offset + (2*U0+2-1-j)*(j)/2 + (l-j-1)];
				}
			}
			sum += m->data[i][m->width - 1];
			sum %= 2;
			mm.data[i][RSD_K - RSD_H - F * U - U0] = sum;
		}

		if (dbg) {
			cout << "dbg final mm:" << endl;
			mm.printIt();
		}


		auto gmm = justGaussLeftUp(mm,mm.width-1);
		if (dbg) {
			cout << " dbg after gauss:" << endl;
			gmm.printIt();
		}

		// check
		bool check_flag = true;
		if (gmm.height == 1) {
			cout << "Final step: gauss fail!!!" << endl;
		}
		for (int i = 0; i < RSD_K - RSD_H - F * U - U0; i++) {
			if (err2.data[i + RSD_N - RSD_K] != gmm.data[i][gmm.width - 1])
			{
				check_flag = false;
				break;
			}
		}
		if (check_flag && cheat_ok) {
			cout << "Finally Verified Okay!(knowing the answer)" << endl;

			for (int i = 0; i < RSD_N; i++)partAns.data[i] = 0;

			for (int i = 0; i < RSD_K - RSD_H - F * U - U0; i++) {
				partAns.data[i] = gmm.data[i][gmm.width - 1];
			}
			for (int i = 0; i < Q; i++) {
				if(pos[i] < RSD_B-1)
					partAns.data[i * (RSD_B - 1) + RSD_K - RSD_H - F * U - U0 + pos[i]] = 1;
			}
			partAns.data[RSD_K - RSD_H - F * U] = 1; //����
			return 0;
		}
		else {
			//cout << "check fail at here" << endl;
			return -1;
		}

	}
	else {
		for (int i = 0; i < RSD_B; i++) {
			pos[tmpPos] = i;
			if (guessParPos(pos, tmpPos + 1, m) == 0)return 0;;
		}
	}
	return -1;
}

int gauss_rsd_attack2() {

	// generate H e = s

	// init the question
	// N+1 include the syndrome column
	Matrix mm(RSD_N - RSD_K, RSD_N + 1);
	// init it with 0 s.t. when matrix * vec is right
	for (int i = 0; i < RSD_N - RSD_K; i++) {
		mm.data[i][RSD_N] = 0;
	}


	int i = 0;
	// add RSD_H or not

	for (i = 0; i < RSD_N - RSD_K; i++) {
		for (int j = 0; j < RSD_N; j++) {
			if (rand() % 12 < 6)mm.data[i][j] = 0;
			else mm.data[i][j] = 1;

		}
	}



	// N+1 include 1 as constant
	
	int err_idx[1000];
	for (int i = 0; i < RSD_N; i++)err.data[i] = 0;
	for (int i = 0; i < RSD_H; i++) {
		int r = rand() % RSD_B;
		err_idx[i] = r;
		// no F U trick
		//if (i < F) {
		//	r = rand() % (RSD_B - U) + U;
		//}

		err.data[i * RSD_B + r] = 1;
	}
	err.data[RSD_N] = 1;
	cout << "err vec:" << endl; err.printIt(RSD_B);
	if (dbg) { cout << "err vec:" << endl; err.printIt(RSD_B); }

	Vec synd = matrixMult(mm, err);
	if (dbg) { cout << "synd vec:" << endl; synd.printIt(); }


	for (int i = 0; i < RSD_N - RSD_K; i++) {
		mm.data[i][RSD_N] = synd.data[i];
	}

	if (dbg) { cout << "orig mat:" << endl;  mm.printIt(); }

	// N -> N-H
	Matrix mm1(RSD_N - RSD_K, RSD_N + 1 - RSD_H);
	for (int i = 0; i < RSD_N - RSD_K; i++) {
		// construct the i-th row of mm1

		// copy s[i]
		mm1.data[i][RSD_N - RSD_H] = mm.data[i][RSD_N];
		for (int j = 0; j < RSD_H; j++) {
			char flag = mm.data[i][j * RSD_B + RSD_B - 1];
			for (int k = 0; k < RSD_B - 1; k++) {
				if (flag)mm1.data[i][j * (RSD_B - 1) + k] = mm.data[i][j * (RSD_B)+k] ^ 1;
				else mm1.data[i][j * (RSD_B - 1) + k] = mm.data[i][j * (RSD_B)+k];
			}
			if (flag)mm1.data[i][RSD_N - RSD_H] ^= 1;
		}
	}
	if (dbg) {
		cout << "mm1" << endl;
		mm1.printIt();
	}


	bool error_free = false;
	int removedColumns[RSD_H][RSD_B - 1];
	
	
	
	do {
		error_free = true;
		for (int i = 0; i < RSD_H; i++) {
			for (int j = 0; j < RSD_B - 1; j++) {
				removedColumns[i][j] = 1;
			}
		}

		for (int i = 0; i < F; i++) {
			int sum = 0;
			while (sum < U) {
				int l = rand() % (RSD_B - 1);
				if (removedColumns[i][l] == 1) {
					removedColumns[i][l] = 0;
					sum++;
					if (err.data[i * RSD_B + l] == 1) {
						error_free = false;
					}
				}

			}
		}
		// just print if we guess fu correct
		if (error_free) {
			cout << "(cheat) This time fu guess error free" << endl;
		}
		else {
			cout << "Guess fu again" << endl;
			//continue;
		}

		// copy the remaining index
		int remainNums[F][RSD_B - 1 - U];
		for (int i = 0; i < F; i++) {
			int idx = 0;
			for (int j = 0; j < RSD_B - 1; j++) {
				if (removedColumns[i][j] == 1) {
					remainNums[i][idx] = j;
					idx++;
				}
			}
		}
		//dbg
		if (dbg) {
			cout << "dbg fu choose:" << endl;
			for (int i = 0; i < F; i++) {
				for (int j = 0; j < RSD_B - 1; j++) {
					cout << removedColumns[i][j];
				}cout << " ";
			}cout << endl;
		}

		// refresh f,u
		Matrix mm2(RSD_N - RSD_K, RSD_N - RSD_H - F * U + 1);

		for (int i = 0; i < RSD_N - RSD_H + 1 - F * U; i++)err2.data[i] = 0;

		for (int i = 0; i < F; i++) {
			for (int j = 0; j < RSD_B - U - 1; j++) {
				//special
				err2.data[i * (RSD_B - U - 1) + j] = err.data[i * RSD_B + remainNums[i][j]];
			}
		}
		for (int i = 0; i < RSD_H - F; i++) {
			for (int j = 0; j < RSD_B - 1; j++) {
				err2.data[F * (RSD_B - U - 1) + i * (RSD_B - 1) + j] = err.data[(F + i) * RSD_B + j];
			}
		}

		err2.data[RSD_N - RSD_H - F * U] = err.data[RSD_N];
		if (dbg) {
			cout << "err2" << endl;
			err2.printIt();
		}

		for (int i = 0; i < RSD_N - RSD_K; i++) {
			for (int j = 0; j < F; j++) {
				for (int k = 0; k < RSD_B - U - 1; k++) {
					mm2.data[i][j * (RSD_B - U - 1) + k] = mm1.data[i][j * (RSD_B - 1) + remainNums[j][k]];
				}
			}
			int offset = F * (RSD_B - U - 1);
			for (int j = 0; j < RSD_H - F; j++) {
				for (int k = 0; k < RSD_B - 1; k++) {
					mm2.data[i][offset + j * (RSD_B - 1) + k] = mm1.data[i][F * (RSD_B - 1) + j * (RSD_B - 1) + k];
				}
			}
			mm2.data[i][RSD_N - RSD_H - F * U] = mm1.data[i][RSD_N - RSD_H];
		}

		if (dbg) {
			cout << "mm2" << endl;
			mm2.printIt();
		}


		// Gauss 1
		Matrix gaussed_mm2 = justGauss2(mm2);
		if (gaussed_mm2.height == 1) {
			//return -1;
			// give fu another choice
			continue;
		}
		if (dbg) {
			cout << "gaussed mat:" << endl; gaussed_mm2.printIt();
		}

		if (important_print) {
			cout << "first time gaussed over" << endl;
		}

		const int m0 = (RSD_H - F - Q) * (RSD_B - 1) * (RSD_B - 2) / 2 + F * (RSD_B - 1 - U) * (RSD_B - 2 - U) / 2;
		const int v0 = RSD_K - RSD_H - F * U - U0;
		const int var_num = v0;
		const int par_num = U0 + 1; // treat 1 as special parameter


		//cout << "n0=" << RSD_K - RSD_H - F * U << endl;
		//cout << "u0=" << U0 << endl;

		ParMat parMat(m0, var_num, par_num); //dbg
		int tmpIdx = 0;

		char bgVecInit[v0];
		char edVecInit[v0];
		char bgParInit[par_num];	// all zero is ok
		char edParInit[par_num];

		for (int iii = 0; iii < par_num; iii++)bgParInit[iii] = edParInit[iii] = 0;

		for (int block = 0; block < F; block++) {
			for (int i = 0; i < RSD_B - U - 1; i++) {
				int bgIdx = block * (RSD_B - U - 1) + i;
				// get the bg Vector
				char* bgVec;
				if (bgIdx < RSD_N - RSD_K) bgVec = &gaussed_mm2.data[bgIdx][RSD_N - RSD_K];
				else {
					for (int iii = 0; iii < v0; iii++)bgVecInit[iii] = 0;
					bgVecInit[bgIdx - RSD_N + RSD_K] = 1;
					bgVec = &bgVecInit[0];
				}
				char* bgPar;
				if (bgIdx < RSD_N - RSD_K) bgPar = &gaussed_mm2.data[bgIdx][RSD_N - RSD_K + var_num];
				else {
					bgPar = &bgParInit[0];
				}

				for (int j = i + 1; j < RSD_B - U - 1; j++) {
					// get the ed Vector
					int edIdx = block * (RSD_B - U - 1) + j;
					char* edVec;
					if (edIdx < RSD_N - RSD_K) edVec = &gaussed_mm2.data[edIdx][RSD_N - RSD_K];
					else {
						for (int iii = 0; iii < v0; iii++)edVecInit[iii] = 0;
						edVecInit[edIdx - RSD_N + RSD_K] = 1;
						edVec = &edVecInit[0];
					}
					char* edPar;
					if (edIdx < RSD_N - RSD_K) edPar = &gaussed_mm2.data[edIdx][RSD_N - RSD_K + var_num];
					else {
						edPar = &edParInit[0];
					}


					parMat.addRow(tmpIdx++, bgVec, bgPar, edVec, edPar);
				}
			}
		}

		for (int block = F; block < RSD_H - Q; block++) {
			for (int i = 0; i < RSD_B - 1; i++) {
				int bgIdx = F * (RSD_B - U - 1) + (block - F) * (RSD_B - 1) + i;
				// get the bg Vector
				char* bgVec;
				if (bgIdx < RSD_N - RSD_K) bgVec = &gaussed_mm2.data[bgIdx][RSD_N - RSD_K];
				else {
					for (int iii = 0; iii < v0; iii++)bgVecInit[iii] = 0;
					bgVecInit[bgIdx - RSD_N + RSD_K] = 1;
					bgVec = &bgVecInit[0];
				}
				char* bgPar;
				if (bgIdx < RSD_N - RSD_K) bgPar = &gaussed_mm2.data[bgIdx][RSD_N - RSD_K + var_num];
				else {
					bgPar = &bgParInit[0];
				}

				for (int j = i + 1; j < RSD_B - 1; j++) {
					// get the ed Vector
					int edIdx = F * (RSD_B - U - 1) + (block - F) * (RSD_B - 1) + j;
					char* edVec;
					if (edIdx < RSD_N - RSD_K) edVec = &gaussed_mm2.data[edIdx][RSD_N - RSD_K];
					else {
						for (int iii = 0; iii < v0; iii++)edVecInit[iii] = 0;
						edVecInit[edIdx - RSD_N + RSD_K] = 1;
						edVec = &edVecInit[0];
					}
					char* edPar;
					if (edIdx < RSD_N - RSD_K) edPar = &gaussed_mm2.data[edIdx][RSD_N - RSD_K + var_num];
					else {
						edPar = &edParInit[0];
					}
					if (dbg) {
						cout << "dbg:" << tmpIdx << endl;
						for (int iii = 0; iii < var_num; iii++) {
							cout << (int)bgVec[iii] << " ";
						}cout << endl;
						for (int iii = 0; iii < var_num; iii++) {
							cout << (int)edVec[iii] << " ";
						}cout << endl;
					}
					parMat.addRow(tmpIdx++, bgVec, bgPar, edVec, edPar);
				}
			}
		}


		parMat.gaussPrepare();
		auto gm = parMat.gaussedParM;

		if(dbg)cout << "!!!!!! ready to guess the real parameter" << endl;
		// guess the real parameter
		int pos[Q];
		int guessFigureok = guessParPos(pos, 0, gm);
		if (guessFigureok == 0) {
			if (vrfy_dbg) {
				cout << "dbg partAns" << endl;
				partAns.printIt();
			}


			Vec fuFreeAns(RSD_N);
			for (int i = 0; i < RSD_N - RSD_K; i++) {
				int sum = 0;
				for (int j = 0; j < RSD_K - RSD_H - F * U + 1; j++) {
					sum += gaussed_mm2.data[i][RSD_N - RSD_K + j] * partAns.data[j];
				}
				sum %= 2;
				fuFreeAns.data[i] = sum;
			}
			for (int i = RSD_N - RSD_K; i < RSD_N - RSD_H - F * U + 1; i++) {
				fuFreeAns.data[i] = partAns.data[i - RSD_N + RSD_K];
			}
			if (vrfy_dbg) {
				cout << "dbg fuFreeAns:" << endl;
				fuFreeAns.printIt();
			}

			Vec hFreeAns(RSD_N);
			int srcIt = 0;

			for (int i = 0; i < F; i++) {
				for (int j = 0; j < RSD_B - 1; j++) {
					hFreeAns.data[i * (RSD_B - 1) + j] = 0;
				}

				for (int j = 0; j < RSD_B - 1 - U; j++) {
					hFreeAns.data[i * (RSD_B - 1) + remainNums[i][j]] = fuFreeAns.data[srcIt++];
				}
			}
			for (int i = F * (RSD_B - 1); i < RSD_N - RSD_H + 1; i++)hFreeAns.data[i] = fuFreeAns.data[srcIt++];

			if (vrfy_dbg) {
				cout << "dbg hFreeAns:" << endl;
				hFreeAns.printIt();
			}

			Vec finalAns(RSD_N + 1);
			for (int i = 0; i < RSD_H; i++) {
				int sum = 0;
				for (int j = 0; j < RSD_B - 1; j++) {
					finalAns.data[i * RSD_B + j] = hFreeAns.data[i * (RSD_B - 1) + j];
					sum += hFreeAns.data[i * (RSD_B - 1) + j];
				}
				sum = 1 - sum;
				sum %= 2;
				finalAns.data[i * RSD_B + RSD_B - 1] = sum;
			}

			if (vrfy_dbg) {
				cout << "dbg finalAns:" << endl;
				finalAns.printIt();
				cout << "dbg err:" << endl;
				err.printIt();
			}

			bool finalAnsOk = true;
			for (int i = 0; i < RSD_N; i++) {
				if (err.data[i] != finalAns.data[i]) {
					finalAnsOk = false;
					break;
				}
			}
			if (finalAnsOk)cout << "FINAL ANS OKAY!!!" << endl;
			else cout << "THIS ANS NOT OK???" << endl;

			if (finalAnsOk)break;
		}
	} while (true);


	
	
	
	return 0;


}

int parMatrixTest() {
	int varNum = 5;
	int parNum = 6;
	char v1[5] = {0,0,1,1,0};
	char v2[5] = { 1,0,0,1,0 };
	char p1[6] = { 0,0,1,1,0,1 };
	char p2[6] = { 1,1,1,0,0,1 };
	auto pm = ParMat(2, 5, 6);
	//pm.printRow(0);
	pm.addRow(0, v1, p1, v2, p2);
	pm.printRow(0);
	pm.addRow(1, v2, p1, v1, p2);
	pm.printRow(1);
	pm.gaussPrepare();
	return 0;
}

int full_attack() {
	int seed = (unsigned)time(NULL);
	srand(seed);
	cout << "bg" << endl;
	cout << "n=" << RSD_N << " k=" << RSD_K << " h=" << RSD_H << endl;
	cout << "f=" << F << " u=" << U << " g=" << Q << endl;

	cout << "m=" << F * (RSD_B - 1 - U) * (RSD_B - 2 - U) / 2 + (RSD_H - F - Q) * (RSD_B - 1) * (RSD_B - 2) / 2 << endl;
	cout << "v=" << RSD_K - RSD_H - F * U - U0 << endl;
	int v = RSD_K - RSD_H - F * U - U0;
	cout << "(v^2+v)/2=" << (v * v + v) / 2 << endl;


	int i = 0;
	while (gauss_rsd_attack2() == -1) {
		cout << "trying:" << i++ << endl;
	}
	cout << "at " << i << " iter finally ok" << endl;
	return 0;
}

int main() {
	//parMatrixTest();
	full_attack();
	system("pause");
	return 0;
}
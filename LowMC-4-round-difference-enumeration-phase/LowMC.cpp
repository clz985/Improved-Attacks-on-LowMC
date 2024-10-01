#include "LowMC.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

#include <random>
#include <ctime>
using namespace std;

double C(int n, int m)
{
	if (m < n - m) m = n - m;
	double ans = 1;
	for (int i = m + 1; i <= n; i++) ans *= i;
	for (int i = 1; i <= n - m; i++) ans /= i;
	return ans;
}
LowMC::LowMC() {
	L = new M[R];
	IL = new M[R];
	KF = new M[R + 1];
	CONS.resize(R);
	for (int i = 0; i < R; i++) {
		CONS[i].resize(N);
	}
	//load round constantC:\Users\anny\Desktop\lowMC_algebraicMITM-main1\lowMC_algebraicMITM-main\constant.txt
	ifstream con("constantFull.txt");
	for (int i = 0; i < R; i++) {
		bool a = 0;
		for (int j = 0; j < N; j++) {
			con >> a;
			CONS[i][j] = a;
		}
	}
	con.close();
	//load linear layers
	initializePars("linearFull.txt", 0, L);

	//load the inverse of the linear layers
	initializePars("inverseFull.txt", 0, IL);
	//load the key schedule function
	initializePars("keyFull.txt", 1, KF);

	//bool checknon[N] = {
	//	0,1,1,1,1,0,0,1,1,1,1,0,1,1,1,0,1,1,0,1,0,1,0,0,1,1,0,1,1,0,0,1,0,1,0,1,0,1,1,0,0,0,1,1,1,1,0,1,1,1,0,1,1,1,1,1,0,0,1,1,1,1,1,1,1,0,0,1,0,1,0,0,0,1,0,1,1,1,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,1,0,0,0,1,0,1,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,1,1,0,1,0,0,
	//};

	/*M id;
	id.r = N;
	id.c = N;
	clearMatrix(id);
	matrixMul(L[102], IL[102], id);
	outputM(id);*/

	bool out[N];
	//	matrixMul(IL[102], checknon, out);
		//bool gamma[N] = {
		//	0,1,0,1,0,1,0,0,1,0,1,0,1,1,0,1,0,1,1,0,0,1,1,1,0,0,1,1,0,0,0,1,0,1,0,0,0,1,1,0,0,1,1,0,1,0,1,1,0,0,1,0,0,1,1,1,1,1,1,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,0,0,1,1,0,1,0,0,0,1,0,0,1,0,0,0,0,1,0,1,0,1,0,1,1,1,1,1,0,1,0,0,0,0,0,1,0,0,1,1,1,0,0,1,0,0,1,1,0,1,0,1,1,0,
		//};
		/*
		for (int i = 0; i < N; i++) {
			if (gamma[i] != out[i]) {
				cout << "wrong1" << endl;
			}
		}
		*/
}

LowMC::~LowMC() {
	delete[]L;
	delete[]IL;
	delete[]KF;
}

void LowMC::initializePars(string filename, bool isKey, M* mat) {
	ifstream in(filename);
	int length = R;
	if (isKey)
		length++;
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				in >> mat[i].ma[j][k];
				//cout << mat[i].ma[j][k];
			}
		}
	}
	in.close();
}

/*****************************
*the encryption function
******************************/
void LowMC::encryptionDetails(bool p0[], bool k[], bool c[], int rounds, vector<vector<bool> >& LOut, vector<vector<bool> >& SOut) {
	M rk;
	rk.r = R + 1;
	rk.c = N;

	bool* p1;
	p1 = new bool[N];
	for (int i = 0; i < N; i++) {
		p1[i] = p0[i];
	}

	for (int i = 0; i < R + 1; i++) {
		matrixMul(KF[i], k, rk.ma[i]);
	}
	//whitening key
	for (int i = 0; i < N; i++) {
		p1[i] = p1[i] ^ rk.ma[0][i];
	}

	int t[3];
	bool* pt = new bool[N];
	for (int i = 0; i < rounds; i++) {
		//s-box (only works for the first 3m bits)
		for (int j = 0; j < N / 3; j++) {
			t[0] = p1[3 * j] ^ (p1[1 + 3 * j] & p1[2 + 3 * j]);
			t[1] = p1[3 * j] ^ p1[1 + 3 * j] ^ (p1[3 * j] & p1[2 + 3 * j]);
			t[2] = p1[3 * j] ^ p1[1 + 3 * j] ^ p1[2 + 3 * j] ^ (p1[3 * j] & p1[1 + 3 * j]);
			p1[3 * j] = t[0];
			p1[3 * j + 1] = t[1];
			p1[3 * j + 2] = t[2];
		}
		//the value after SBox
		for (int j = 0; j < N; j++) {
			SOut[i][j] = p1[j];
		}

		//linear matrix
		matrixMul(L[i], p1, pt);

		//constant addition
		for (int j = 0; j < N; j++) {
			p1[j] = pt[j] ^ CONS[i][j];
		}

		//key addition
		for (int j = 0; j < N; j++) {
			p1[j] = p1[j] ^ rk.ma[i + 1][j];
		}

		for (int j = 0; j < N; j++) {
			LOut[i][j] = p1[j];
		}
	}
	for (int i = 0; i < N; i++) {
		c[i] = p1[i];
	}
	delete[]pt;
	delete[]p1;
}

/*****************************
*
*used for attacks
*
*****************************/




/********************************
*constructEquation
*********************************/





void LowMC::constructEquation2(vector<vector<bool> > SOutDiff, vector<vector<bool> > LOutDiff, int Eqnumber, int g1, int g2)
{
	int Snum = 2 * g2;
	int ccol = 2 * (N / 3 - g1) + 3 * g2 + 1;
	int col = ccol * (ccol - 1) * (ccol - 2) / 6 + ccol * (ccol - 1) / 2 + ccol + 1;
	M E1;
	E1.r = N;
	E1.c = ccol;
	clearMatrix(E1);
	M E1E;
	E1E.r = 3 * (N / 3 - g1);
	E1E.c = ccol;
	clearMatrix(E1E);

	bool input1[N];
	for (int i = 0; i < E1E.r; i++)
	{
		input1[i] = LOutDiff[Eqnumber - 1][i + 3 * g1];
	}
	constructExpressions(input1, E1E, ccol);
	
	for (int i = 0; i < 3 * g1; i++)
	{
		for (int j = 0; j < E1.c; j++)
		{
			E1.ma[i][j] = 0;
		}
		E1.ma[i][E1.c - 1] = SOutDiff[Eqnumber][i];

	}
	for (int i = 0; i < E1E.r; i++)
	{
		for (int j = 0; j < E1E.c; j++)
		{
			E1.ma[3 * g1 + i][j] = E1E.ma[i][j];
		}
	}
	//expression of the input in the second round 
	matrixMul(L[Eqnumber], E1);
	
	M E2;
	E2.r = N;
	E2.c = ccol;
	clearMatrix(E2);
	for (int i = 0; i < g2; i++)
	{
		E2.ma[3 * i][3 * i + 2 * (N / 3 - g1)] = 1;
		E2.ma[3 * i + 1][3 * i + 1 + 2 * (N / 3 - g1)] = 1;
		E2.ma[3 * i + 2][3 * i + 2 + 2 * (N / 3 - g1)] = 1;
	}
	vector<vector<bool> > cubic;
	vector<vector<bool> > Cubic;
	cubic.resize(100);
	Cubic.resize(100);
	for (int i = 0; i < 100; i++) {
		cubic[i].resize(col);
		Cubic[i].resize(col);
	}
	for (int i = 0; i < cubic.size(); i++) {
		for (int j = 0; j < cubic[i].size(); j++) {
			cubic[i][j] = 0;
		}
	}
	for (int i = 0; i < Cubic.size(); i++) {
		for (int j = 0; j < Cubic[i].size(); j++) {
			Cubic[i][j] = 0;
		}
	}


	for (int k = 0; k < g2; k++)
	{

		

		addLinear1(cubic, 2 * k, E1, 3 * k, col);


		addLinear1(cubic, 2 * k, E1, 3 * k + 1, col);


		addLinear1(cubic, 2 * k, E1, 3 * k + 2, col);


		addQuadratic1(cubic, 2 * k, E1, 3 * k, E1, 3 * k + 1, ccol - 1, col);


		addQuadratic1(cubic, 2 * k, E1, 3 * k, E1, 3 * k + 2, ccol - 1, col);


		addQuadratic1(cubic, 2 * k, E1, 3 * k + 1, E1, 3 * k + 2, ccol - 1, col);


		addCubic1(cubic, 2 * k, E1, 3 * k, E1, 3 * k + 1, E1, 3 * k + 2, ccol - 1, col);



		addLinear1(cubic, 2 * k + 1, E1, 3 * k, col);

		addLinear1(cubic, 2 * k + 1, E1, 3 * k + 1, col);

		addLinear1(cubic, 2 * k + 1, E1, 3 * k + 2, col);

		addQuadratic1(cubic, 2 * k + 1, E1, 3 * k, E1, 3 * k + 1, ccol - 1, col);

		addQuadratic1(cubic, 2 * k + 1, E1, 3 * k, E1, 3 * k + 2, ccol - 1, col);

		addQuadratic1(cubic, 2 * k + 1, E1, 3 * k + 1, E1, 3 * k + 2, ccol - 1, col);

		addCubic1(cubic, 2 * k + 1, E1, 3 * k, E1, 3 * k + 1, E1, 3 * k + 2, ccol - 1, col);
	}

	cout << "Cubicleft:" << endl;
	for (int i = 0; i < Snum; i++) {
		for (int j = 0; j < cubic[i].size(); j++) {
			cout << cubic[i][j];
		}
		cout << endl;
	}
	cout << endl;




	

	for (int k = 0; k < g2; k++)
	{



		addLinear1(Cubic, 2 * k, E2, 3 * k, col);
		addLinear1(Cubic, 2 * k, E2, 3 * k + 1, col);
		addLinear1(Cubic, 2 * k, E2, 3 * k + 2, col);
		addQuadratic1(Cubic, 2 * k, E2, 3 * k, E2, 3 * k + 1, ccol - 1, col);
		addQuadratic1(Cubic, 2 * k, E2, 3 * k, E2, 3 * k + 2, ccol - 1, col);
		addQuadratic1(Cubic, 2 * k, E2, 3 * k + 1, E2, 3 * k + 2, ccol - 1, col);
		addCubic1(Cubic, 2 * k, E2, 3 * k, E2, 3 * k + 1, E2, 3 * k + 2, ccol - 1, col);

		addQuadratic1(Cubic, 2 * k + 1, E1, 3 * k, E2, 3 * k, ccol - 1, col);
		addQuadratic1(Cubic, 2 * k + 1, E1, 3 * k + 1, E2, 3 * k + 1, ccol - 1, col);
		addQuadratic1(Cubic, 2 * k + 1, E1, 3 * k + 2, E2, 3 * k + 2, ccol - 1, col);
	}

	cout << "Cubicright:" << endl;
	for (int i = 0; i < Snum; i++) {
		for (int j = 0; j < Cubic[i].size(); j++) {
			cout << Cubic[i][j];
		}
		cout << endl;
	}
	cout << endl;

	for (int i = 0; i < cubic.size(); i++) {
		cubic[i].clear();
	}
	cubic.clear();
	for (int i = 0; i < Cubic.size(); i++) {
		Cubic[i].clear();
	}
	Cubic.clear();

}

void LowMC::constructExpressions(bool diff[], M& eq, int col) {
	int varNum = 0;

	bool* con = new bool[N];//the constant part in the expression
	for (int i = 0; i < N; i++) {
		con[i] = 0;
	}
	clearMatrix(eq);

	int inactive = 0;
	//assign variables based on diff
	for (int i = 0; i < eq.r / 3; i++) {//deal with each s-box
		int a = diff[3 * i];
		int b = diff[3 * i + 1];
		int c = diff[3 * i + 2];
		a = a + b * 2 + c * 4;
		if (a == 0) {//(0,0,0)
			//do nothing
			inactive++;
		}
		else if (a == 4) {//(0,0,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			varNum += 2;//two variables are added
		}
		else if (a == 2) {//(0,1,0)
			eq.ma[3 * i][varNum] = 1;
			con[3 * i + 1] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			varNum += 2;
		}
		else if (a == 1) {//(1,0,0)
			con[3 * i] = 1;
			eq.ma[3 * i + 1][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			varNum += 2;
		}
		else if (a == 6) {//(0,1,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			varNum += 2;
		}
		else if (a == 3) {//(1,1,0)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum] = 1;
			con[3 * i + 1] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			varNum += 2;
		}
		else if (a == 5) {//(1,0,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			con[3 * i + 2] = 1;
			varNum += 2;
		}
		else {//(1,1,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			varNum += 2;
		}
	}
	//cout << "varNum: " << varNum << endl;
	//cout << "inactive:" << inactive << endl;
	//change eq
	//eq.c = varNum + 1;
	for (int i = 0; i < N; i++) {
		eq.ma[i][eq.c - 1] = con[i];
	}


	delete[]con;
}

void LowMC::addLinear1(vector<vector<bool>>& qeq, int row, M& lin, int index, int col) {
	for (int i = 0; i < lin.c - 1; i++) {//variable part
		qeq[row][i] = qeq[row][i] ^ lin.ma[index][i];
	}
	qeq[row][col - 1] = qeq[row][col - 1] ^ lin.ma[index][lin.c - 1];//constant part
}

void LowMC::addQuadratic1(vector<vector<bool>>& qeq, int row, M& fac0, int in0, M& fac1, int in1, int total, int col) {
	//(f0+a) & (f1+b) = f0f1 + b*f0 + a*f1 + ab
	if (fac0.ma[in0][fac0.c - 1]) {//a=1, add f1
		for (int i = 0; i < fac1.c - 1; i++) {//variable part
			qeq[row][i] = qeq[row][i] ^ fac1.ma[in1][i];
		}
	}

	if (fac1.ma[in1][fac1.c - 1]) {
		for (int i = 0; i < fac0.c - 1; i++) {//variable part
			qeq[row][i] = qeq[row][i] ^ fac0.ma[in0][i];
		}
	}

	qeq[row][col - 1] = qeq[row][col - 1] ^ (fac0.ma[in0][fac0.c - 1] & fac1.ma[in1][fac1.c - 1]);

	//add f0*f1
	for (int i = 0; i < fac0.c - 1; i++) {
		for (int j = 0; j < fac1.c - 1; j++) {
			if (fac0.ma[in0][i] == 1 && fac1.ma[in1][j] == 1) {
				//the term x_i*x_j
				int i0 = i, j0 = j;//(j0>=i0 is necessary)
				if (i0 > j0) {
					swap(i0, j0);
				}
				int index = (total + (total - (j0 - i0) + 1)) * (j0 - i0) / 2 + i0;
				if (index >= col) {
					cout << "i,j:" << i << " , " << j << endl;
					cout << "wrong8" << endl;
				}
				if (index < 0) {
					cout << "minus" << endl;
				}
				if (index < col && index >= 0)
				{
					qeq[row][index] = qeq[row][index] ^ 1;
				}
			}
		}
	}
}

void LowMC::addCubic1(vector<vector<bool> >& qeq, int row, M& fac0, int in0, M& fac1, int in1, M& fac2, int in2, int total, int col)
{//(f0+a)&(f1+b)&(f2+c)=f0f1f2+c*f0f1+b*f0f2+a*f1f2+ab*f2+ac*f1+bc*f0+abc
	if (fac0.ma[in0][fac0.c - 1])
	{//a=1,add f1f2
		for (int i = 0; i < fac1.c - 1; i++) {
			for (int j = 0; j < fac2.c - 1; j++) {
				if (fac1.ma[in1][i] == 1 && fac2.ma[in2][j] == 1) {
					//the term x_i*x_j
					int i0 = i, j0 = j;//(j0>=i0 is necessary)
					if (i0 > j0) {
						swap(i0, j0);
					}
					int index = (total + (total - (j0 - i0) + 1)) * (j0 - i0) / 2 + i0;
					if (index >= total * (total - 1) / 2) {
						cout << "i,j:" << i << " , " << j << endl;
						cout << "wrong1" << endl;
					}
					if (index < 0) {
						cout << "minus" << endl;
					}
					qeq[row][index] = qeq[row][index] ^ 1;
				}
			}
		}

	}

	if (fac1.ma[in1][fac1.c - 1]) {//b=1,add f0f2,f2
		for (int i = 0; i < fac0.c - 1; i++) {
			for (int j = 0; j < fac2.c - 1; j++) {
				if (fac0.ma[in0][i] == 1 && fac2.ma[in2][j] == 1) {
					//the term x_i*x_j
					int i0 = i, j0 = j;//(j0>=i0 is necessary)
					if (i0 > j0) {
						swap(i0, j0);
					}
					int index = (total + (total - (j0 - i0) + 1)) * (j0 - i0) / 2 + i0;

					if (index >= total * (total - 1) / 2) {
						cout << "i,j:" << i << " , " << j << endl;
						cout << "wrong2" << endl;
					}
					if (index < 0) {
						cout << "minus" << endl;
					}
					qeq[row][index] = qeq[row][index] ^ 1;
				}
			}
		}
		for (int i = 0; i < fac2.c - 1; i++) {//variable part
			qeq[row][i] = qeq[row][i] ^ fac2.ma[in2][i];
		}
	}

	if (fac2.ma[in2][fac1.c - 1])
	{//c=1,add f0f1,f1,f0
		for (int i = 0; i < fac0.c - 1; i++) {
			for (int j = 0; j < fac1.c - 1; j++) {
				if (fac0.ma[in0][i] == 1 && fac1.ma[in1][j] == 1) {
					//the term x_i*x_j
					int i0 = i, j0 = j;//(j0>=i0 is necessary)
					if (i0 > j0) {
						swap(i0, j0);
					}
					int index = (total + (total - (j0 - i0) + 1)) * (j0 - i0) / 2 + i0;

					if (index >= total * (total - 1) / 2) {
						cout << "i,j:" << i << " , " << j << endl;
						cout << "wrong3" << endl;
					}
					if (index < 0) {
						cout << "minus" << endl;
					}
					qeq[row][index] = qeq[row][index] ^ 1;
				}
			}
		}


		for (int i = 0; i < fac1.c - 1; i++) {//variable part
			qeq[row][i] = qeq[row][i] ^ fac1.ma[in1][i];
		}

		for (int i = 0; i < fac0.c - 1; i++) {//variable part
			qeq[row][i] = qeq[row][i] ^ fac0.ma[in0][i];
		}

	}

	qeq[row][col - 1] = qeq[row][col - 1] ^ (fac0.ma[in0][fac0.c - 1] & fac1.ma[in1][fac1.c - 1] & fac2.ma[in2][fac2.c - 1]);
	//add f0*f1*f2

	for (int i = 0; i < fac0.c - 1; i++) {
		for (int j = 0; j < fac1.c - 1; j++) {
			for (int k = 0; k < fac2.c - 1; k++) {
				if (fac0.ma[in0][i] == 1 && fac1.ma[in1][j] == 1) {
					if (fac2.ma[in2][k] == 1)
					{
						//the term x_i*x_j*x_k
						int i0 = i, j0 = j, k0 = k;//(j0>=i0 is necessary)




						int m;
						if (i0 > j0) {
							m = i0;
							i0 = j0;
							j0 = m;
						}
						if (j0 > k0) {
							m = j0;
							j0 = k0;
							k0 = m;
						}
						if (j0 < i0) {
							m = i0;
							i0 = j0;
							j0 = m;
						}






						int r = j0 - i0;

						int f = k0 - j0;



						if (r == 0)
						{
							if (f > 0)
							{
								int index = (total + (total - (k0 - i0) + 1)) * (k0 - i0) / 2 + i0;
								qeq[row][index] = qeq[row][index] ^ 1;

							}
							else
							{
								int index = i0;
								qeq[row][index] = qeq[row][index] ^ 1;

							}
						}
						else
						{
							if (f == 0)
							{
								int index = (total + (total - (k0 - i0) + 1)) * (k0 - i0) / 2 + i0;
								qeq[row][index] = qeq[row][index] ^ 1;

							}
							else {
								int w = 0;
								for (int i = 0; i < i0; i++)
								{
									w += C(total - 1 - i, 2);
								}
								int index = total + C(total, 2) + w + (total - i0 - 2 + (total - i0 - 1 - (k0 - j0) + 1)) * (k0 - j0 - 1) / 2 + (j0 - i0 - 1);

								if (index >= col) {
									cout << "i,j,k:" << i << " , " << j << " , " << k << endl;
									cout << "wrong4" << endl;
								}
								if (index < 0) {
									cout << "minus1" << endl;
								}
								qeq[row][index] = qeq[row][index] ^ 1;

							}

						}


					}



				}
			}
		}
	}






}



/********************************
*the optimized key-recovery phase
*********************************/


void LowMC::keyRecovery(bool key[], bool c[], int r, vector<vector<bool> >& LOut, vector<vector<bool> >& SOut, vector<vector<bool> >& sOut, vector<vector<bool> >& lOut, bool y[]) {
	//Round function = S -> L
	int total = 3;//consider the last 3 rounds
	int S = 118;//consider the last 118 S-boxes
	int index = 0;
	vector<vector<bool> > intermediateSOut(total);
	for (int i = 0; i < total; i++) {
		intermediateSOut[i].clear();
		for (int j = 0; j < N; j++)
			intermediateSOut[i].push_back(0);
	}
	bool tmp[N];
	for (int i = 0; i < N; i++)
		tmp[i] = c[i];

	M eqSys;
	eqSys.r = 354;//118*3 (at most 243 linear equations)
	eqSys.c = 199;//129+3*23+1
	clearMatrix(eqSys);

	vector<M> expLOut(total + 1);
	for (int i = 0; i < expLOut.size(); i++) {
		expLOut[i].r = N, expLOut[i].c = 199;
		clearMatrix(expLOut[i]);
	}
	for (int i = 0; i < N; i++)
		expLOut[0].ma[i][expLOut[0].c - 1] = c[i];//initialize the constant part

	vector<M> expSOut(total + 1);
	for (int i = 0; i < expSOut.size(); i++) {
		expSOut[i].r = N, expSOut[i].c = 199;
		clearMatrix(expSOut[i]);
	}

	int cnt = 0;
	int activenumer = 0;
	eqSys.r = 0;//at first, there is no useful equation.
	//test the correctness of expressions
	bool testVec[199], resVec[199];
	for (int i = 0; i < 199; i++)
		testVec[i] = 0;
	for (int i = 0; i < N; i++)
		testVec[i] = key[i];
	testVec[198] = 1;
	/////////////////////////////////////

	int freeVars = 0;//the number of free variables after gaussian elimination
	//store the expressions of the intermediate variables
	M extraEq;
	extraEq.r = 69;//23*3
	extraEq.c = 199;
	clearMatrix(extraEq);

	for (int i = 0; i < total; i++) {
		index = r - 1 - i;
		//first, inverse constant addition
		for (int j = 0; j < N; j++)
			expLOut[i].ma[j][expLOut[i].c - 1] ^= CONS[index][j];
		//add round key KF[index+1]
		for (int jr = 0; jr < N; jr++) {
			for (int jc = 0; jc < N; jc++)
				expLOut[i].ma[jr][jc] ^= KF[index + 1].ma[jr][jc];
		}
		//inverse L
		matrixMul(IL[index], expLOut[i], expSOut[i]);
		///////////check the correctness of expSOut[i]
		/*
		cout << "check the correctness of expSOut[i]" << endl;
		matrixMul(expSOut[i], testVec, resVec);
		for (int j = 0; j < 129; j++)
			resVec[j] ^= SOut[index][j];
		outputArray(resVec, 129);
		*/

		///////////////////////
		//inverse S (use sOut[index],lOut[index-1])
		//check whether the output difference is zero
		//cout << "out=" << endl;
		for (int Sbox = 0; Sbox < 43; Sbox++)
		{
			//cout << "Sbox=" << Sbox << endl;
			int out = sOut[index][0 + 3 * Sbox] + 2 * sOut[index][1 + 3 * Sbox] + 4 * sOut[index][2 + 3 * Sbox];
			int in = lOut[index - 1][0 + 3 * Sbox] + 2 * lOut[index - 1][1 + 3 * Sbox] + 4 * lOut[index - 1][2 + 3 * Sbox];

			//cout << hex << in << " " << out << endl;

			if (i < 2)
			{
				if (out == 0) {// intermeidate variables
					expLOut[i + 1].ma[0 + 3 * Sbox][N + 3 * cnt] = 1;
					expLOut[i + 1].ma[1 + 3 * Sbox][N + 3 * cnt + 1] = 1;
					expLOut[i + 1].ma[2 + 3 * Sbox][N + 3 * cnt + 2] = 1;
					//update testVintroduceector
					testVec[N + 3 * cnt] = LOut[index - 1][0 + 3 * Sbox];
					testVec[N + 3 * cnt + 1] = LOut[index - 1][1 + 3 * Sbox];
					testVec[N + 3 * cnt + 2] = LOut[index - 1][2 + 3 * Sbox];
					///////////////////
					//store the expressions
					for (int jr = 0; jr < 3; jr++) {
						for (int jc = 0; jc < expSOut[i].c; jc++)
							extraEq.ma[3 * cnt + jr][jc] = expSOut[i].ma[jr + 3 * Sbox][jc];
					}
					cnt++;
				}

				else {//extract equations and freely linearize the sbox
					dynamicallyUpdateExpression(expLOut[i + 1], expSOut[i], in, out, Sbox);
					extractEquations(eqSys, expSOut[i], in, out, Sbox);
				}

			}
			else
			{
				if (out == 0) continue;
				else
				{
					extractEquations(eqSys, expSOut[i], in, out, Sbox);
				}
			}
			/*
			////check the correctness of the extract equations
			cout << "check the correctness of the extract equations" << endl;
			if (eqSys.r > 0) {
				matrixMul(eqSys, testVec, resVec);
				outputArray(resVec, eqSys.r);
				cout << "i=" << i << endl;
				if (resVec[eqSys.r - 1] != 0 || resVec[eqSys.r - 2] != 0) {
					system("pause");
				}
			}
			///////the correctness is verified
			*/
		}
		/*
		if (i < 2)
		{
			//cout << "expLOut[i + 1]" << endl;
			//outputM(expLOut[i + 1]);
			//outputArray(testVec, 199);
			
			////////////check the correctness of expLOut[i+1]
			matrixMul(expLOut[i + 1], testVec, resVec);
			for (int j = 0; j < 129; j++)
				resVec[j] ^= LOut[index - 1][j];
			cout << "check the correctness of expLOut[i+1]" << endl;
			outputArray(resVec, 129);


			if (resVec[0] != 0 || resVec[1] != 0 || resVec[2] != 0) {
				cout << SOut[index][0] << " " << SOut[index][1] << " " << SOut[index][2] << endl;
				cout << LOut[index - 1][0] << " " << LOut[index - 1][1] << " " << LOut[index - 1][2] << endl;
				//outputM(expSOut[i]);
				//cout << i<<"diff:" << in << " " << out << endl;
				//outputM(expLOut[i + 1]);
				system("pause");
			}
		}
		//////////////////////////////////////////////////the correctness is verified
		*/


		//the condition to exit
			if (eqSys.r >= N + 3 * cnt) {//directly solve the system
				cout << "linear:" << eqSys.r << " \t quadratic:" << 0 << endl;
				break;
			}

		}
		cout << "cnt=" << cnt << endl;
		cout << "eqSys.r=" << eqSys.r << endl;

		//outputM(extraEq);
		freeVars = (N + 3 * cnt) - eqSys.r;
		if (freeVars > 0) {
			cout << "eqSys.r=" << eqSys.r << endl;
			cout << "freeVars:" << freeVars << endl;
			cout << "row:" << cnt * 3 << endl;
		}


		bool yy[5];
		for (int i = 0; i < 5; i++)
		{
			yy[i] = y[i];
		}

		solveKey(eqSys, extraEq, cnt * 3, testVec, key, yy);

		//system("pause");
		//cout << "the size of the constructed equation system:" << eqSys.r << endl;

		//gauss(eqSys, cnt * 3 + 128);
		//outputM(eqSys);
		//system("pause");
		//outputM(expLOut[81]);
		//system("pause");
	
}

	void LowMC::solveKey(M& eqSys, M& extraEq, int extraSize, bool testVec[], bool key[], bool yy[]) {
		/*cout << "test" << endl;

		for (int i = 0; i < N + extraSize; i++) {
			cout << testVec[i];
		}
		cout << endl;
		cout << "eqSys=" << endl;
		outputM(eqSys);
		cout << endl;*/
		gauss(eqSys, extraSize + 129);
		//outputM(eqSys);
		cout << "further simplified:" << endl;
		simplify(eqSys, extraSize + 129);
		//outputM(eqSys);
		if (eqSys.r >= N + extraSize)
		{//directly solve the system
			cout << "directly";
			return;
		}
		if (eqSys.r < N + extraSize)
		{

			/*
			cout << "check:";
			bool res[199];
			matrixMul(eqSys, testVec, res);
			
			for (int i = 0; i < eqSys.r; i++)
				cout << res[i];
			cout << endl;
			*/
			/*cout << endl;
			cout << "eqSys=" << endl;
			outputM(eqSys);
			cout << endl;*/
			//mark the free variables
			int totalVars = extraSize + 129;
			vector<bool> freeVar(totalVars);
			vector<int> index(totalVars);
			for (int i = 0; i < totalVars; i++)
				index[i] = 0;
			for (int i = 0; i < totalVars; i++)
				freeVar[i] = 1;
			for (int i = 0; i < eqSys.r; i++) {
				for (int j = i; j < eqSys.c - 1; j++) {
					if (eqSys.ma[i][j] == 1) {
						freeVar[j] = 0;
						index[j] = i;//record the expression for the variable u_j, i.e. the i-th row
						break;
					}
				}
			}
			/*vector<bool> zz(totalVars);
			cout << "zz=" << endl;
			for (int i = 0; i < totalVars; i++)
			{zz[i] = 0;

				if (freeVar[i]) {
					zz[i] = 1;
				}
				cout << zz[i] ;
			}
			cout << endl;*/
			int freeVarCnt = 0;
			for (int i = 0; i < freeVar.size(); i++) {
				if (freeVar[i])
					freeVarCnt++;
			}
			cout << "the number of free variables:" << freeVarCnt << " \t #terms:" << (freeVarCnt + (freeVarCnt) * (freeVarCnt - 1) / 2) << endl;

			//update extraEq
			for (int i = 0; i < extraSize; i++) {
				for (int j = 0; j < extraEq.c - 1; j++) {
					if (extraEq.ma[i][j] == 1 && freeVar[j] == 0) {
						extraEq.ma[i][j] = 0;//update it
						for (int k = j + 1; k < extraEq.c; k++) {
							extraEq.ma[i][k] = extraEq.ma[i][k] ^ eqSys.ma[index[j]][k];
						}
					}
				}
			}

			//outputM(extraEq);
			/*matrixMul(extraEq, testVec, res);
			bool t[3];
			for (int j = 0; j < extraSize / 3; j++) {
				t[0] = testVec[129+3 * j] ^ (testVec[129+1 + 3 * j] & testVec[129+2 + 3 * j]);
				t[1] = testVec[129+3 * j] ^ testVec[129+1 + 3 * j] ^ (testVec[129+3 * j] & testVec[129+2 + 3 * j]);
				t[2] = testVec[129+3 * j] ^ testVec[129+1 + 3 * j] ^ testVec[129+2 + 3 * j] ^ (testVec[129+3 * j] & testVec[129+1 + 3 * j]);
				cout << t[0] << t[1] << t[2];
			}
			cout << endl;
			outputArray(res, extraSize);*/
			///the correctness is verified

			//rename variables
			vector<int> mapIndex(totalVars), inverseMap(totalVars);
			int rvar = 0;
			for (int i = 0; i < totalVars; i++) {
				if (freeVar[i]) {
					mapIndex[i] = rvar;
					inverseMap[rvar] = i;
					rvar++;
				}
			}

			M left, right;
			left.r = extraSize, right.r = extraSize;
			left.c = freeVarCnt + 1, right.c = freeVarCnt + 1;
			clearMatrix(left), clearMatrix(right);

			bool* checkVec = new bool[freeVarCnt + 1];
			for (int i = 0; i < freeVarCnt; i++) {
				checkVec[i] = testVec[inverseMap[i]];

			}
			checkVec[freeVarCnt] = 1;
			/*
			cout << "checkVec" << endl;
			for (int i = 0; i < freeVarCnt + 1; i++)
			{
				cout << checkVec[i];
			}
			cout << endl;
			*/
			for (int i = 0; i < extraSize; i++) {
				if (freeVar[129 + i]) {//it is a free variables
					left.ma[i][mapIndex[129 + i]] = 1;
					left.ma[i][left.c - 1] = 0;
				}
				else {//it is not a free variable, get the corresponding expression
					for (int j = 129 + i + 1; j < eqSys.c - 1; j++) {
						if (eqSys.ma[index[129 + i]][j] == 1) {
							left.ma[i][mapIndex[j]] = 1;
							//cout << j << " " << mapIndex[j] << endl;
						}
					}
					left.ma[i][left.c - 1] = eqSys.ma[index[129 + i]][eqSys.c - 1];
				}
			}
			//cout << "left:" << endl;
			//outputM(left);//test is passed, left is correctly constructed

			//construct right
			for (int i = 0; i < extraSize; i++) {
				for (int j = 0; j < extraEq.c - 1; j++) {
					if (extraEq.ma[i][j] == 1) {
						right.ma[i][mapIndex[j]] = 1;
					}
				}
				right.ma[i][right.c - 1] = extraEq.ma[i][extraEq.c - 1];
			}
			//cout << "right:" << endl;
			//outputM(right);//test is passed, left is correctly constructed

			//construct the quadratic equation system with left and right
			int quTerms = freeVarCnt + (freeVarCnt * (freeVarCnt - 1) / 2);
			bool* quVec = new bool[quTerms + 1];
			for (int i = 0; i < freeVarCnt; i++) {
				for (int j = i; j < freeVarCnt; j++) {
					int k = (freeVarCnt + (freeVarCnt - (j - i) + 1)) * (j - i) / 2 + i;
					if (k >= quTerms) {
						cout << "wrong5" << endl;
						system("pause");
					}
					quVec[k] = checkVec[i] & checkVec[j];
					//cout << checkVec[i] << checkVec[j] << quVec[k] << " ";
				}
			}
			//cout << endl;
			quVec[quTerms] = 1;
			/*cout << "queVec" << endl;
			for (int i = 0; i < quTerms; i++) {
				cout << quVec[i];
			}*/
			//cout << endl;
			bool quVecRes[300];
			bool* quu = new bool[7];
			for (int i = 0; i < 3; i++) {
				for (int j = i + 1; j < 3; j++) {
					if (i < freeVarCnt && j < freeVarCnt)
					{
						int k = (3 + (3 - (j - i) + 1)) * (j - i) / 2 + i;

						quu[k - 3] = checkVec[i] & checkVec[j];
						//cout << checkVec[i] << checkVec[j] << quVec[k] << " ";
					}
				}
			}
			for (int i = 0; i < 3; i++)
			{
				if (i < freeVarCnt)
				{
					quu[i + 3] = checkVec[i];
				}
			}
			quu[6] = 1;
			//construct quadratic equations
			M quadratic;
			quadratic.r = 1 * 14;
			quadratic.c = quTerms + 1;
			clearMatrix(quadratic);

			//cout << quadratic.r << " " << quadratic.c << endl;
			//cout << "quadratic:" << endl;
			//outputM(quadratic);
			//cout << "updated:" << endl;

			int cnt = 1;
			for (int i = 0; i < cnt; i++) {
				//left.ma[3*i], ma[3*i+1], ma[3*i+2]
				//right.ma[3*i], ma[3*i+1], ma[3*i+2]
				int j = 3 * i;
				int k = 14 * i;

				addLinear(quadratic, k, left, j);
				addLinear(quadratic, k, right, j);
				addQuadratic(quadratic, k, left, j + 1, left, j + 2, freeVarCnt);
				//matrixMul(quadratic, quVec, quVecRes);
				//cout << "index:" << k - 14 * i << ": ";
				//outputArray(quVecRes, 37);
				k++;

				addLinear(quadratic, k, left, j);
				addLinear(quadratic, k, left, j + 1);
				addLinear(quadratic, k, right, j + 1);
				addQuadratic(quadratic, k, left, j, left, j + 2, freeVarCnt);
				//matrixMul(quadratic, quVec, quVecRes);
				//cout << "index:" << k - 14 * i << ": ";
				//outputArray(quVecRes, 37);
				k++;

				addLinear(quadratic, k, left, j);
				addLinear(quadratic, k, left, j + 1);
				addLinear(quadratic, k, left, j + 2);
				addLinear(quadratic, k, right, j + 2);
				addQuadratic(quadratic, k, left, j, left, j + 1, freeVarCnt);
				//matrixMul(quadratic, quVec, quVecRes);
				//cout << "index:" << k - 14 * i << ": ";
				//outputArray(quVecRes, 37);
				k++;

				addLinear(quadratic, k, right, j);
				addLinear(quadratic, k, right, j + 1);
				addLinear(quadratic, k, left, j);
				addQuadratic(quadratic, k, right, j + 1, right, j + 2, freeVarCnt);
				//matrixMul(quadratic, quVec, quVecRes);
				//cout << "index:" << k - 14 * i << ": ";
				//outputArray(quVecRes, 37);
				k++;

				addLinear(quadratic, k, right, j + 1);
				addLinear(quadratic, k, left, j + 1);
				addQuadratic(quadratic, k, right, j, right, j + 2, freeVarCnt);
				//matrixMul(quadratic, quVec, quVecRes);
				//cout << "index:" << k - 14 * i << ": ";
				//outputArray(quVecRes, 37);
				k++;

				addLinear(quadratic, k, right, j);
				addLinear(quadratic, k, right, j + 1);
				addLinear(quadratic, k, right, j + 2);
				addLinear(quadratic, k, left, j + 2);
				addQuadratic(quadratic, k, right, j, right, j + 1, freeVarCnt);
				//matrixMul(quadratic, quVec, quVecRes);
				//cout << "index:" << k - 14 * i << ": ";
				//outputArray(quVecRes, 37);
				k++;

				addQuadratic(quadratic, k, right, j, left, j + 1, freeVarCnt);
				addQuadratic(quadratic, k, left, j, left, j + 1, freeVarCnt);
				addQuadratic(quadratic, k, left, j + 1, left, j + 2, freeVarCnt);
				//matrixMul(quadratic, quVec, quVecRes);
				//cout << "index:" << k - 14 * i << ": ";
				//outputArray(quVecRes, 37);
				k++;

				addQuadratic(quadratic, k, right, j, left, j + 2, freeVarCnt);
				addQuadratic(quadratic, k, left, j, left, j + 2, freeVarCnt);
				addQuadratic(quadratic, k, left, j + 1, left, j + 2, freeVarCnt);
				//matrixMul(quadratic, quVec, quVecRes);
				//cout << "index:" << k - 14 * i << ": ";
				//outputArray(quVecRes, 37);
				k++;

				addQuadratic(quadratic, k, right, j + 1, left, j, freeVarCnt);
				addLinear(quadratic, k, left, j);
				addQuadratic(quadratic, k, left, j, left, j + 1, freeVarCnt);
				addQuadratic(quadratic, k, left, j, left, j + 2, freeVarCnt);
				//matrixMul(quadratic, quVec, quVecRes);
				//cout << "index:" << k - 14 * i << ": ";
				//outputArray(quVecRes, 37);
				k++;

				addQuadratic(quadratic, k, right, j + 1, left, j + 2, freeVarCnt);
				addQuadratic(quadratic, k, left, j + 1, left, j + 2, freeVarCnt);
				//matrixMul(quadratic, quVec, quVecRes);
				//cout << "index:" << k - 14 * i << ": ";
				//outputArray(quVecRes, 37);
				k++;

				//system("pause");
				//eq11,12,13,14
				addQuadratic(quadratic, k, right, j + 2, left, j, freeVarCnt);
				addLinear(quadratic, k, left, j);
				addQuadratic(quadratic, k, left, j, left, j + 2, freeVarCnt);
				//matrixMul(quadratic, quVec, quVecRes);
				//cout << "index:" << k - 14 * i << ": ";
				//outputArray(quVecRes, 37);
				k++;

				addQuadratic(quadratic, k, right, j + 2, left, j + 1, freeVarCnt);
				addLinear(quadratic, k, left, j + 1);
				addQuadratic(quadratic, k, left, j + 1, left, j + 2, freeVarCnt);
				//matrixMul(quadratic, quVec, quVecRes);
				//cout << "index:" << k - 14 * i << ": ";
				//outputArray(quVecRes, 37);
				k++;

				addQuadratic(quadratic, k, right, j, left, j, freeVarCnt);
				addLinear(quadratic, k, left, j);
				addQuadratic(quadratic, k, right, j + 1, left, j + 1, freeVarCnt);
				addQuadratic(quadratic, k, left, j, left, j + 1, freeVarCnt);
				addLinear(quadratic, k, left, j + 1);
				//matrixMul(quadratic, quVec, quVecRes);
				//cout << "index:" << k - 14 * i << ": ";
				//outputArray(quVecRes, 37);
				k++;

				addQuadratic(quadratic, k, right, j + 1, left, j + 1, freeVarCnt);
				addQuadratic(quadratic, k, left, j, left, j + 1, freeVarCnt);
				addLinear(quadratic, k, left, j + 1);
				addQuadratic(quadratic, k, right, j + 2, left, j + 2, freeVarCnt);
				addQuadratic(quadratic, k, left, j, left, j + 2, freeVarCnt);
				addQuadratic(quadratic, k, left, j + 1, left, j + 2, freeVarCnt);
				addLinear(quadratic, k, left, j + 2);
				//matrixMul(quadratic, quVec, quVecRes);
				//cout << "index:" << k - 14 * i << ": ";
				//outputArray(quVecRes, 37);
				k++;
			}

			//matrixMul(quadratic, quVec, quVecRes);
			//for (int i = 0; i < quadratic.c; i++) {
			//	if (quVecRes[i] != 0) {
			//		cout << "wrong6" << endl;
			//		system("pause");
			//	}
			//}
			//cout << "quadratic=" << endl;
			//outputM(quadratic);
			//cout << endl;
			M Updatequadratic;
			Updatequadratic.r = 1 * 14;
			Updatequadratic.c = quTerms + 1;
			//x0*x1--x0*x2

			// x0*x1--x1*x2

			for (int i = 0; i < Updatequadratic.r; i++)
			{
				for (int j = 0; j <= 1; j++)
				{
					Updatequadratic.ma[i][j] = quadratic.ma[i][8 + j];
				}
			}
			//x0*x2
			for (int i = 0; i < Updatequadratic.r; i++)
			{
				for (int j = 0; j <= 0; j++)
				{
					Updatequadratic.ma[i][j + 2] = quadratic.ma[i][15 + j];
				}
			}



			//l(y)*x0

			for (int i = 0; i < Updatequadratic.r; i++)
			{
				int w0 = 15;//15+0
				for (int j = 6; j > 1; j--)
				{


					w0 = j + w0;
					Updatequadratic.ma[i][3 + 6 - j] = quadratic.ma[i][w0];


				}
			}
			//l(y)*x1

			for (int i = 0; i < Updatequadratic.r; i++)
			{
				int w1 = 9;//8+1
				for (int j = 7; j > 2; j--)
				{


					w1 = j + w1;
					Updatequadratic.ma[i][7 - j + 8] = quadratic.ma[i][w1];


				}
			}
			//l(y)*x2

			for (int i = 0; i < Updatequadratic.r; i++)
			{
				int w2 = 2;//2
				for (int j = 8; j > 3; j--)
				{


					w2 = j + w2;
					Updatequadratic.ma[i][8 - j + 13] = quadratic.ma[i][w2];


				}
			}

			//x0--x2
			for (int i = 0; i < Updatequadratic.r; i++)
			{
				for (int j = 0; j <= 2; j++)
				{
					Updatequadratic.ma[i][j + 18] = quadratic.ma[i][j];
				}
			}
			//x3x4--x3x7
			//x3x4--x6x7
			for (int i = 0; i < Updatequadratic.r; i++)
			{
				for (int j = 0; j <= 3; j++)
				{
					Updatequadratic.ma[i][j + 21] = quadratic.ma[i][8 + 3 + j];
				}
			}
			//x3x5--x5x7
			for (int i = 0; i < Updatequadratic.r; i++)
			{
				for (int j = 0; j <= 2; j++)
				{
					Updatequadratic.ma[i][j + 25] = quadratic.ma[i][15 + 3 + j];
				}
			}
			//x3x6--x4x7
			for (int i = 0; i < Updatequadratic.r; i++)
			{
				for (int j = 0; j <= 1; j++)
				{
					Updatequadratic.ma[i][j + 28] = quadratic.ma[i][21 + 3 + j];
				}
			}
			//x3x7
			for (int i = 0; i < Updatequadratic.r; i++)
			{
				for (int j = 0; j <= 0; j++)
				{
					Updatequadratic.ma[i][j + 30] = quadratic.ma[i][26 + 3 + j];
				}
			}

			//x3--x7
			for (int i = 0; i < Updatequadratic.r; i++)
			{
				for (int j = 0; j <= 4; j++)
				{
					Updatequadratic.ma[i][j + 31] = quadratic.ma[i][j + 3];
				}
			}
			for (int i = 0; i < Updatequadratic.r; i++)
			{

				Updatequadratic.ma[i][Updatequadratic.c - 1] = quadratic.ma[i][quadratic.c - 1];

			}
			//cout << "Updatequadratic=" << endl;
			//outputM(Updatequadratic);
			//cout << endl;
			//Adjust the order to fit the crossbred algorithm
			M colExQuadratic;
			colExQuadratic.r = Updatequadratic.r;
			colExQuadratic.c = Updatequadratic.c;
			clearMatrix(colExQuadratic);
			/*
			for (int i = 0; i < 5; i++)
			{
				cout << yy[i];
			}
			cout << endl;
			*/
			for (int i = 0; i < Updatequadratic.r; i++) {
				for (int j = 0; j < 3; j++) {
					colExQuadratic.ma[i][j] = Updatequadratic.ma[i][j];
				}
				colExQuadratic.ma[i][colExQuadratic.c - 1] = Updatequadratic.ma[i][Updatequadratic.c - 1];
			}
			for (int i = 0; i < Updatequadratic.r; i++) {
				for (int j = 0; j < 3; j++) {
					colExQuadratic.ma[i][18 + j] = Updatequadratic.ma[i][18 + j];
				}
			}

			for (int k = 0; k < Updatequadratic.r; k++) {
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 5; j++) {

						if (yy[j] == 1)
						{
							if (Updatequadratic.ma[k][3 + i * 5 + j] == 1) {
								colExQuadratic.ma[k][3 + i * 5 + j] = 1;
							}

						}
						else {
							colExQuadratic.ma[k][3 + i * 5 + j] = 0;
						}
					}
				}
			}

			for (int k = 0; k < Updatequadratic.r; k++) {
				for (int i = 0; i < 5; i++) {
					if (yy[i] == 1 && Updatequadratic.ma[k][31 + i] == 1)
					{
						colExQuadratic.ma[k][31 + i] = 1;
					}
				}

			}

			for (int k = 0; k < Updatequadratic.r; k++) {
				for (int i = 0; i < 5; i++) {
					for (int j = i + 1; j < 5; j++) {
						int s = 0;
						s = (5 + (5 - (j - i) + 1)) * (j - i) / 2 + i - 3 - 5 + 21;
						if (yy[i] == 1 && yy[j] == 1)
						{
							if (Updatequadratic.ma[k][s] == 1) {
								colExQuadratic.ma[k][s] = 1;
							}
						}

					}
				}
			}
			//cout << "colExQuadratic=" << endl;
			//outputM(colExQuadratic);
			//cout << endl;

			M col;
			col.r = colExQuadratic.r;
			col.c = 3 + colExQuadratic.r;
			clearMatrix(col);
			for (int i = 0; i < colExQuadratic.r; i++) {
				for (int j = 0; j < 3; j++) {
					col.ma[i][j] = colExQuadratic.ma[i][j];
				}
				col.ma[i][3 + i] = 1;
			}
			//cout << "Quadratic coefficient matrix=" << endl;
			//outputM(col);
			//cout << endl;
			gauss(col, col.c);
			//simplify(col, col.c);
			
			//cout << "Quadratic coefficient matrix=" << endl;
			//outputM(col);
			//cout << endl;
			
			M pp;
			pp.r = col.r;
			pp.c = col.r;
			clearMatrix(pp);
			for (int i = 0; i < col.r; i++) {

				for (int j = 0; j < col.r; j++) {
					pp.ma[i][j] = col.ma[i][j + 3];
				}
			}
			//cout << "pp=" << endl;
			//outputM(pp);
			//cout << endl;
			M colExQuadratic1;
			colExQuadratic1.r = 1 * 14;
			colExQuadratic1.c = 7;
			//l(y)*x0
			for (int i = 0; i < colExQuadratic1.r; i++) {
				int z = 0;
				for (int j = 0; j < 5; j++) {
					z = z ^ colExQuadratic.ma[i][3 + j];
				}
				colExQuadratic1.ma[i][3] = z ^ colExQuadratic.ma[i][18];
			}
			//l(y)*x1
			for (int i = 0; i < colExQuadratic1.r; i++) {


				int z = 0;
				for (int j = 0; j < 5; j++) {
					z = z ^ colExQuadratic.ma[i][8 + j];

				}
				colExQuadratic1.ma[i][4] = z ^ colExQuadratic.ma[i][19];

			}
			//l(y)*x2
			for (int i = 0; i < colExQuadratic1.r; i++) {


				int z = 0;
				for (int j = 0; j < 5; j++) {
					z = z ^ colExQuadratic.ma[i][13 + j];

				}
				colExQuadratic1.ma[i][5] = z ^ colExQuadratic.ma[i][20];

			}



			for (int i = 0; i < colExQuadratic1.r; i++) {


				int z = 0;
				for (int j = 21; j < colExQuadratic.c; j++) {
					z = z ^ colExQuadratic.ma[i][j];
				}


				colExQuadratic1.ma[i][6] = z;

			}
			for (int i = 0; i < colExQuadratic1.r; i++) {



				for (int j = 0; j < 3; j++) {
					colExQuadratic1.ma[i][j] = colExQuadratic.ma[i][j];
				}
			}

			//test have pasted
			/*

			matrixMul(colExQuadratic1, quu, quVecRes);
			cout << "quVecRes=" << endl;
			for (int i = 0; i < 18; i++)
			{
				cout << quVecRes[i] ;
			}
			cout << endl;*/
			//cout << "colExQuadratic1=" << endl;
			//outputM(colExQuadratic1);
			//cout << endl;

			M qq;
			qq.r = colExQuadratic1.r;
			qq.c = 4;
			for (int i = 0; i < colExQuadratic1.r; i++) {
				for (int j = 0; j < 4; j++) {
					qq.ma[i][j] = colExQuadratic1.ma[i][j + 3];
				}
			}
			matrixMul(pp, qq);
			//cout << "qq=" << endl;
			//outputM(qq);
			//cout << endl;

			//cout << "solve last 11 variables" << endl;
			M solve;
			solve.r = 14;
			solve.c = 4;
			for (int i = 0; i < 14; i++) {
				for (int j = 0; j < 4; j++) {
					solve.ma[i][j] = qq.ma[i][j];
				}
			}


			//cout << "solve=" << endl;
			//outputM(solve);
			//cout << endl;
			gauss(solve, solve.c);

			//cout << "solve=" << endl;
			//outputM(solve);
			//cout << endl;
			vector<vector<bool> > sol;
			int solNum = 0;

			storeSolutions(sol, solve, solNum);
			cout << "solNum=" << solNum << endl;
			/*

			cout << "checkVec=" << endl;
			for (int i = 0; i < 15; i++) {
				cout << checkVec[i];
			}
			cout << endl;
			*/



			for (int i = 0; i < solNum; i++) {
				bool findKey = true;
				for (int j = 0; j < 3; j++) {
					//cout << sol[i][j];
					if (sol[i][j] != checkVec[j]) {

						findKey = false;
					}

				}
				cout << endl;
				if (findKey)
				{			
							cout << "sol" << i << " : correct" << endl;
					
				}
					
				else
					cout << "sol" << i << " :wrong" << endl;
				//outputArray(checkVec, freeVarCnt);
			}

			cout << "solve end" << endl;
			cout << endl;

			delete[]checkVec;
			delete[]quVec;
			//system("pause");
		}
	
}

void LowMC::addLinear(M& qeq, int row, M& lin, int index) {
	for (int i = 0; i < lin.c - 1; i++) {//variable part
		qeq.ma[row][i] ^= lin.ma[index][i];
	}
	qeq.ma[row][qeq.c - 1] ^= lin.ma[index][lin.c - 1];//constant part
}

void LowMC::addQuadratic(M& qeq, int row, M& fac0, int in0, M& fac1, int in1, int total) {
	//(f0+a) & (f1+b) = f0f1 + b*f0 + a*f1 + ab
	if (fac0.ma[in0][fac0.c - 1]) {//a=1, add f1
		for (int i = 0; i < fac1.c - 1; i++) {//variable part
			qeq.ma[row][i] ^= fac1.ma[in1][i];
		}
	}

	if (fac1.ma[in1][fac1.c - 1]) {
		for (int i = 0; i < fac0.c - 1; i++) {//variable part
			qeq.ma[row][i] ^= fac0.ma[in0][i];
		}
	}

	qeq.ma[row][qeq.c - 1] = qeq.ma[row][qeq.c - 1] ^ (fac0.ma[in0][fac0.c - 1] & fac1.ma[in1][fac1.c - 1]);

	//add f0*f1
	for (int i = 0; i < fac0.c - 1; i++) {
		for (int j = 0; j < fac1.c - 1; j++) {
			if (fac0.ma[in0][i] == 1 && fac1.ma[in1][j] == 1) {
				//the term x_i*x_j
				int i0 = i, j0 = j;//(j0>=i0 is necessary)
				if (i0 > j0) {
					int t = i0;
					i0 = j0;
					j0 = t;
				}
				int index = (total + (total - (j0 - i0) + 1)) * (j0 - i0) / 2 + i0;
				if (index >= qeq.c) {
					cout << "i,j:" << i << " , " << j << endl;
					cout << "wrong8" << endl;
				}
				if (index < 0) {
					cout << "minus" << endl;
				}
				if (index < qeq.c && index >= 0)
				{
					qeq.ma[row][index] ^= 1;
				}

			}
		}
	}
}

void LowMC::dynamicallyUpdateExpression(M& expLOut, M& expSOut, int in, int out, int Sbox) {
	bool c[3] = { 0,0,0 };
	vector<vector<int> > add(3);
	for (int i = 0; i < 3; i++)
		add[i].clear();

	//input diff : x0x1x2: 100
	if (in == 1 && out == 1) {
		//x0=z0+1, x1=1, x2=1, 
		c[0] = 1, c[1] = 1, c[2] = 1;
		add[0].push_back(0);
	}
	else if (in == 1 && out == 5) {
		//x0=z0,x1=0,x2=1
		c[0] = 0, c[1] = 0, c[2] = 1;
		add[0].push_back(0);
	}
	else if (in == 1 && out == 3) {
		//x0=z0,x1=1,x2=0
		c[0] = 0, c[1] = 1, c[2] = 0;
		add[0].push_back(0);
	}
	else if (in == 1 && out == 7) {
		//x0=z0,x1=0,x2=0,
		c[0] = 0, c[1] = 0, c[2] = 0;
		add[0].push_back(0);
	}

	//input diff : x0x1x2: 010
	else if (in == 2 && out == 2) {
		//x0=1,x1=z1+1,x2=0
		c[0] = 1, c[1] = 1, c[2] = 0;
		add[1].push_back(1);
	}
	else if (in == 2 && out == 6) {
		//x0=0,x1=z1,x2=0
		c[0] = 0, c[1] = 0, c[2] = 0;
		add[1].push_back(1);
	}
	else if (in == 2 && out == 3) {
		//x0=1,x1=z1,x2=1
		c[0] = 1, c[1] = 0, c[2] = 1;
		add[1].push_back(1);
	}
	else if (in == 2 && out == 7) {
		//x0=0,x1=z1,x2=1
		c[0] = 0, c[1] = 0, c[2] = 1;
		add[1].push_back(1);
	}

	//input diff: x0x1x2: 110
	else if (in == 3 && out == 2) {
		//x0=x1+1, x1=z1, x2=1
		c[0] = 1, c[1] = 0, c[2] = 1;
		add[0].push_back(1), add[1].push_back(1);
	}
	else if (in == 3 && out == 6) {
		//x0=x1, x1=z1, x2=1
		c[0] = 0, c[1] = 0, c[2] = 1;
		add[0].push_back(1), add[1].push_back(1);
	}
	else if (in == 3 && out == 1) {
		//x0=x1+1=z0, x1=z0+1, x2=0
		c[0] = 0, c[1] = 1, c[2] = 0;
		add[0].push_back(0), add[1].push_back(0);
	}
	else if (in == 3 && out == 5) {
		//x0=x1, x1=z0, x2=0
		c[0] = 0, c[1] = 0, c[2] = 0;
		add[0].push_back(0), add[1].push_back(0);
	}

	//input diff: x0x1x2: 001
	else if (in == 4 && out == 4) {
		//x0=0, x1=0, x2=z2
		c[0] = 0, c[1] = 0, c[2] = 0;
		add[2].push_back(2);
	}
	else if (in == 4 && out == 6) {
		//x0=1, x1=0, x2=z2+1
		c[0] = 1, c[1] = 0, c[2] = 1;
		add[2].push_back(2);
	}
	else if (in == 4 && out == 5) {
		//x0=0, x1=1, x2=z2+1
		c[0] = 0, c[1] = 1, c[2] = 1;
		add[2].push_back(2);
	}
	else if (in == 4 && out == 7) {
		//x0=1, x1=1, x2=z2+1
		c[0] = 1, c[1] = 1, c[2] = 1;
		add[2].push_back(2);
	}

	//input diff: x0x1x2: 101
	else if (in == 5 && out == 4) {
		//x0=x2=z2+1, x1=1, x2=z2+1
		c[0] = 1, c[1] = 1, c[2] = 1;
		add[0].push_back(2), add[2].push_back(2);
	}
	else if (in == 5 && out == 1) {
		//x0=x2=z0, x1=0, x2=z0
		c[0] = 0, c[1] = 0, c[2] = 0;
		add[0].push_back(0), add[2].push_back(0);
	}
	else if (in == 5 && out == 6) {
		//x0=x2+1=z1+1, x1=1, x2=z1
		c[0] = 1, c[1] = 1, c[2] = 0;
		add[0].push_back(1), add[2].push_back(1);
	}
	else if (in == 5 && out == 3) {
		//x0=x2+1=z1, x1=0, x2=z1+1
		c[0] = 0, c[1] = 0, c[2] = 1;
		add[0].push_back(1), add[2].push_back(1);
	}

	//input diff: x0x1x2: 011
	else if (in == 6 && out == 4) {
		//x0=1, x1=x2+1, x1=z2 -> x2=z2+1
		c[0] = 1, c[1] = 0, c[2] = 1;
		add[1].push_back(2), add[2].push_back(2);
	}
	else if (in == 6 && out == 2) {
		//x0=0, x1=x2+1, x1=z1 -> x2=z1+1
		c[0] = 0, c[1] = 0, c[2] = 1;
		add[1].push_back(1), add[2].push_back(1);
	}
	else if (in == 6 && out == 5) {
		//x0=1, x1=x2, x1=z0+1 -> x2=z0+1
		c[0] = 1, c[1] = 1, c[2] = 1;
		add[1].push_back(0), add[2].push_back(0);
	}
	else if (in == 6 && out == 3) {
		//x0=0, x1=x2, x1=z0 -> x2=z0
		c[0] = 0, c[1] = 0, c[2] = 0;
		add[1].push_back(0), add[2].push_back(0);
	}

	//input diff: x0x1x2: 111
	else if (in == 7 && out == 4) {
		//x0=z2, x1=x0+1=z2+1, x2=x0+1=z2+1
		c[0] = 0, c[1] = 1, c[2] = 1;
		add[0].push_back(2), add[1].push_back(2), add[2].push_back(2);
	}
	else if (in == 7 && out == 2) {
		//x0=z1, x1=x0=z1, x2=x0=z1
		c[0] = 0, c[1] = 0, c[2] = 0;
		add[0].push_back(1), add[1].push_back(1), add[2].push_back(1);
	}
	else if (in == 7 && out == 1) {
		//x0=z0, x1=x0=z0, x2=x0+1=z0+1
		c[0] = 0, c[1] = 0, c[2] = 1;
		add[0].push_back(0), add[1].push_back(0), add[2].push_back(0);
	}
	else if (in == 7 && out == 7) {
		//x0=z0, x1=x0+1=z0+1, x2=x0=z0
		c[0] = 0, c[1] = 1, c[2] = 0;
		add[0].push_back(0), add[1].push_back(0), add[2].push_back(0);
	}


	for (int i = 0; i < 3; i++) {
		//assign c to the constant part
		expLOut.ma[i + 3 * Sbox][expLOut.c - 1] = c[i];
		//add row (add[i][0]) of expSOut to expLOut
		if (add[i].size() > 0) {
			for (int j = 0; j < expSOut.c; j++)
				expLOut.ma[i + 3 * Sbox][j] ^= expSOut.ma[add[i][0] + 3 * Sbox][j];
		}
	}

}

void LowMC::extractEquations(M& eq, M& expSOut, int in, int out, int Sbox) {
	bool c[2] = { 0,0 };
	vector<vector<int> > add(2);
	add[0].clear(), add[1].clear();

	//input diff : x0x1x2: 100
	if (in == 1 && out == 1) {
		//z1=1, z2=0
		c[0] = 1, c[1] = 0;
		add[0].push_back(1), add[1].push_back(2);
	}
	else if (in == 1 && out == 5) {
		//z1=0, z0+z2=1
		c[0] = 0, c[1] = 1;
		add[0].push_back(1), add[1].push_back(0), add[1].push_back(2);
	}
	else if (in == 1 && out == 3) {
		//z2=1, z0+z1=1
		c[0] = 1, c[1] = 1;
		add[0].push_back(2);
		add[1].push_back(0), add[1].push_back(1);
	}
	else if (in == 1 && out == 7) {
		//z0=z1, z0=z2
		c[0] = 0, c[1] = 0;
		add[0].push_back(0), add[0].push_back(1);
		add[1].push_back(0), add[1].push_back(2);
	}

	//input diff : x0x1x2: 010
	else if (in == 2 && out == 2) {
		//z0=1, z2=1
		c[0] = 1, c[1] = 1;
		add[0].push_back(0);
		add[1].push_back(2);
	}
	else if (in == 2 && out == 6) {
		//z0=0, z1=z2
		c[0] = 0, c[1] = 0;
		add[0].push_back(0);
		add[1].push_back(1), add[1].push_back(2);
	}
	else if (in == 2 && out == 3) {
		//z2=0,z0+z1=1
		c[0] = 0, c[1] = 1;
		add[0].push_back(2);
		add[1].push_back(0), add[1].push_back(1);
	}
	else if (in == 2 && out == 7) {
		//z0+z1=0, z0+z2=1
		c[0] = 0, c[1] = 1;
		add[0].push_back(0), add[0].push_back(1);
		add[1].push_back(0), add[1].push_back(2);
	}

	//input diff: x0x1x2: 110
	else if (in == 3 && out == 2) {
		//z0=1, z2=0
		c[0] = 1, c[1] = 0;
		add[0].push_back(0);
		add[1].push_back(2);
	}
	else if (in == 3 && out == 6) {
		//z0=0, z1+z2=1
		c[0] = 0, c[1] = 1;
		add[0].push_back(0);
		add[1].push_back(1), add[1].push_back(2);
	}
	else if (in == 3 && out == 1) {
		//z1=1, z2=1
		c[0] = 1, c[1] = 1;
		add[0].push_back(1);
		add[1].push_back(2);
	}
	else if (in == 3 && out == 5) {
		//z1=0, z0+z2=0
		c[0] = 0, c[1] = 0;
		add[0].push_back(1);
		add[1].push_back(0), add[1].push_back(2);
	}

	//input diff: x0x1x2: 001
	else if (in == 4 && out == 4) {
		//z0=0, z1=0
		c[0] = 0, c[1] = 0;
		add[0].push_back(0);
		add[1].push_back(1);
	}
	else if (in == 4 && out == 6) {
		//z0=1, z1+z2=0
		c[0] = 1, c[1] = 0;
		add[0].push_back(0);
		add[1].push_back(1), add[1].push_back(2);
	}
	else if (in == 4 && out == 5) {
		//z1=1, z0+z2=1
		c[0] = 1, c[1] = 1;
		add[0].push_back(1);
		add[1].push_back(0), add[1].push_back(2);
	}
	else if (in == 4 && out == 7) {
		//z0+z1=1, z0+z2=0;
		c[0] = 1, c[1] = 0;
		add[0].push_back(0), add[0].push_back(1);
		add[1].push_back(0), add[1].push_back(2);
	}

	//input diff: x0x1x2: 101
	else if (in == 5 && out == 4) {
		//z0=0,z1=1
		c[0] = 0, c[1] = 1;
		add[0].push_back(0);
		add[1].push_back(1);
	}
	else if (in == 5 && out == 1) {
		//z1=0, z2=0
		c[0] = 0, c[1] = 0;
		add[0].push_back(1);
		add[1].push_back(2);
	}
	else if (in == 5 && out == 6) {
		//z0=1, z1+z2=1
		c[0] = 1, c[1] = 1;
		add[0].push_back(0);
		add[1].push_back(1), add[1].push_back(2);
	}
	else if (in == 5 && out == 3) {
		//z2=1, z0+z1=0
		c[0] = 1, c[1] = 0;
		add[0].push_back(2);
		add[1].push_back(0), add[1].push_back(1);
	}

	//input diff: x0x1x2: 011
	else if (in == 6 && out == 4) {
		//z0=1,z1=0
		c[0] = 1, c[1] = 0;
		add[0].push_back(0);
		add[1].push_back(1);
	}
	else if (in == 6 && out == 2) {
		//z0=0, z2=1
		c[0] = 0, c[1] = 1;
		add[0].push_back(0);
		add[1].push_back(2);
	}
	else if (in == 6 && out == 5) {
		//z1=1, z0+z2=0
		c[0] = 1, c[1] = 0;
		add[0].push_back(1);
		add[1].push_back(0), add[1].push_back(2);
	}
	else if (in == 6 && out == 3) {
		//z2=0, z0+z1=0
		c[0] = 0, c[1] = 0;
		add[0].push_back(2);
		add[1].push_back(0), add[1].push_back(1);
	}

	//input diff: x0x1x2: 111
	else if (in == 7 && out == 4) {
		//z0=1,z1=1
		c[0] = 1, c[1] = 1;
		add[0].push_back(0);
		add[1].push_back(1);
	}
	else if (in == 7 && out == 2) {
		//z0=0,z2=0
		c[0] = 0, c[1] = 0;
		add[0].push_back(0);
		add[1].push_back(2);
	}
	else if (in == 7 && out == 1) {
		//z1=0,z2=1
		c[0] = 0, c[1] = 1;
		add[0].push_back(1);
		add[1].push_back(2);
	}
	else if (in == 7 && out == 7) {
		//z0+z1=1, z0+z2=1
		c[0] = 1, c[1] = 1;
		add[0].push_back(0), add[0].push_back(1);
		add[1].push_back(0), add[1].push_back(2);
	}

	//add to eq
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < add[i].size(); j++) {
			for (int k = 0; k < expSOut.c; k++)
				eq.ma[eq.r][k] ^= expSOut.ma[add[i][j] + 3 * Sbox][k];
		}
		eq.ma[eq.r][eq.c - 1] = eq.ma[eq.r][eq.c - 1] ^ c[i];
		eq.r++;//increase the number of eqs
	}
}


/*****************************
*
*algebraic equations for DDT
*
*****************************/

void LowMC::searchAlgebraicEquations() {
	//compute DDT
	int t = 64;//6 bits, x0,x1,x2,x3,x4,x5
	bool x[64][6];
	bool valid[64] = { 0 };
	valid[0] = 1;
	for (int i = 0; i < t; i++) {
		for (int j = 0; j < 6; j++) {
			x[i][j] = (i >> j) & 0x1;
		}
		if ((x[i][0] != 0 || x[i][1] != 0 || x[i][2] != 0) &&
			(x[i][3] != 0 || x[i][4] != 0 || x[i][5] != 0)) {
			bool sum = (x[i][0] & x[i][3]) ^ (x[i][1] & x[i][4]) ^ (x[i][2] & x[i][5]);
			if (sum == 1) {
				valid[i] = 1;
			}
		}
	}

	//compute DDT from the definition
	//x[3]=x[0]+x[1]x[2]
	//x[4]=x[0]+x[1]+x[0]x[2]
	//x[5]=x[0]+x[1]+x[2]+x[0]x[1]
	bool y[6], yp[6];
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			for (int z = 0; z < 8; z++) {
				for (int u = 0; u < 3; u++) {
					y[u] = (z >> u) & 0x1;
					yp[u] = y[u] ^ ((i >> u) & 0x1);
				}
				y[3] = y[0] ^ (y[1] & y[2]);
				y[4] = y[0] ^ y[1] ^ (y[0] & y[2]);
				y[5] = y[0] ^ y[1] ^ y[2] ^ (y[0] & y[1]);

				yp[3] = yp[0] ^ (yp[1] & yp[2]) ^ y[3];
				yp[4] = yp[0] ^ yp[1] ^ (yp[0] & yp[2]) ^ y[4];
				yp[5] = yp[0] ^ yp[1] ^ yp[2] ^ (yp[0] & yp[1]) ^ y[5];

				int out = yp[3] + 2 * yp[4] + 4 * yp[5];
				if (j == out) {
					cout << i << " -> " << yp[3] << " " << yp[4] << " " << yp[5] << endl;
					if (valid[i * 8 + j] == 0) {
						cout << "wrong9" << endl;
					}
					break;
				}
			}
		}
	}

	system("pause");

	//find quadratic equations, we have 6+15+1=18 variables (64 equations)
	M q;
	q.r = 64, q.c = 18;
	clearMatrix(q);
	//outputM(q);

	int c = 0;
	for (int i = 0; i < t; i++)
		q.ma[i][q.c - 1] = 0;//guess the constant

	for (int i = 0; i < t; i++) {
		if (valid[i]) {
			computeQuadraticEq(x[i], q.ma[c]);
			c++;
			//q.ma[i][q.c - 1] = 0;
		}
		else {
			//if ((x[i][0] != 0 || x[i][1] != 0 || x[i][2] != 0) &&
				//(x[i][3] != 0 || x[i][4] != 0 || x[i][5] != 0)) {
			computeQuadraticEq(x[i], q.ma[c]);
			q.ma[c][q.c - 1] ^= 1;
			c++;
			//}
		}
	}
	cout << "c:" << c << endl;
	outputM(q);
	cout << "after gauss" << endl;
	gauss(q, q.c - 1);
	outputM(q);
	cout << "simplify:" << endl;
	simplify(q, q.c - 1);
	outputM(q);

	//find cubic equations, we have 6+15+20+1=42 variables (64 equations)
	q.r = 64, q.c = 42;
	clearMatrix(q);
	c = 0;
	for (int i = 0; i < t; i++)
		q.ma[i][q.c - 1] = 0;//guess the constant
	for (int i = 0; i < t; i++) {
		if (valid[i]) {
			computeCubicEq(x[i], q.ma[c]);
			c++;
			//q.ma[i][q.c - 1] = 0;
		}
		else {
			//if ((x[i][0] != 0 || x[i][1] != 0 || x[i][2] != 0) &&
				//(x[i][3] != 0 || x[i][4] != 0 || x[i][5] != 0)) {
			computeCubicEq(x[i], q.ma[c]);
			q.ma[c][q.c - 1] ^= 1;
			c++;
			//}
		}
	}
	cout << "c:" << c << endl;
	outputM(q);
	cout << "after gauss" << endl;
	gauss(q, q.c - 1);
	outputM(q);
	cout << "simplify:" << endl;
	simplify(q, q.c - 1);
	outputM(q);

	for (int i = 0; i < t; i++) {
		bool c0 = (x[i][0] & x[i][3]) ^ (x[i][1] & x[i][4]) ^ (x[i][2] & x[i][5]);
		c0 = c0 ^ ((x[i][0] ^ 1) & (x[i][1] ^ 1) & (x[i][2] ^ 1));
		c0 = c0 ^ 1;

		bool c1 = ((x[i][0] ^ 1) & (x[i][1] ^ 1) & (x[i][2] ^ 1));
		c1 = c1 ^ ((x[i][3] ^ 1) & (x[i][4] ^ 1) & (x[i][5] ^ 1));

		if (valid[i]) {
			if (c0 != 0 && c1 != 0)
				cout << "wrong valid" << endl;
		}
		else {
			if (c0 == 0 && c1 == 0) {
				cout << x[i][0] << x[i][1] << x[i][2] << " -> " << x[i][3] << x[i][4] << x[i][5] << endl;
				cout << "wrong invalid" << endl;
			}
		}
	}
}

void LowMC::computeQuadraticEq(bool x[], bool row[]) {
	int cnt = 0;
	for (int i = 0; i < 6; i++) {
		for (int j = i + 1; j < 6; j++) {
			row[cnt] = (x[i] & x[j]);
			cnt++;
		}
	}
	for (int i = 0; i < 6; i++) {
		row[cnt] = x[i];
		cnt++;
	}
}

void LowMC::computeCubicEq(bool x[], bool row[]) {
	int cnt = 0;
	for (int i = 0; i < 6; i++) {
		for (int j = i + 1; j < 6; j++) {
			for (int z = j + 1; z < 6; z++) {
				row[cnt] = (x[i] & x[j] & x[z]);
				cnt++;
			}
		}
	}
	for (int i = 0; i < 6; i++) {
		for (int j = i + 1; j < 6; j++) {
			row[cnt] = (x[i] & x[j]);
			cnt++;
		}
	}
	for (int i = 0; i < 6; i++) {
		row[cnt] = x[i];
		cnt++;
	}
}



/*****************************
*
*matrix operations
*
*****************************/

void LowMC::matrixMul(M& m, bool x[], bool y[]) {
	for (int i = 0; i < m.r; i++) {
		y[i] = 0;
		for (int j = 0; j < m.c; j++) {
			y[i] = y[i] ^ (m.ma[i][j] & x[j]);
		}
	}
}

void LowMC::matrixMul(M& m, vector<bool>& x, vector<bool>& y) {
	for (int i = 0; i < m.r; i++) {
		y[i] = 0;
		for (int j = 0; j < m.c; j++) {
			y[i] = y[i] ^ (m.ma[i][j] & x[j]);
		}
	}
}

void LowMC::matrixMul(M& m, bool x[]) {
	bool y[256];
	for (int i = 0; i < m.r; i++) {
		y[i] = 0;
		for (int j = 0; j < m.c; j++) {
			y[i] = y[i] ^ (m.ma[i][j] & x[j]);
		}
	}
	for (int i = 0; i < m.r; i++) {
		x[i] = y[i];
	}
}

void LowMC::matrixMul(M& m1, M& m2, M& m3) {
	//cout << m1.r << " " << m1.c <<" "<<m2.c<< endl;
	for (int i = 0; i < m1.r; i++) {
		for (int j = 0; j < m2.c; j++) {
			m3.ma[i][j] = 0;
			for (int k = 0; k < m1.c; k++) {
				m3.ma[i][j] = m3.ma[i][j] ^ (m1.ma[i][k] & m2.ma[k][j]);
			}
		}
	}
	//outputMatrix(m3);

}

void LowMC::matrixMul(M& m1, M& m2) {
	M m3;
	m3.r = m1.r;
	m3.c = m2.c;
	//cout << m1.r << " " << m1.c <<" "<<m2.c<< endl;
	for (int i = 0; i < m1.r; i++) {
		for (int j = 0; j < m2.c; j++) {
			m3.ma[i][j] = 0;
			for (int k = 0; k < m1.c; k++) {
				m3.ma[i][j] = m3.ma[i][j] ^ (m1.ma[i][k] & m2.ma[k][j]);
			}
		}
	}
	//outputMatrix(m3);
	for (int i = 0; i < m2.r; i++) {
		for (int j = 0; j < m2.c; j++) {
			m2.ma[i][j] = m3.ma[i][j];
		}
	}
}
/*
void LowMC::matrixMul(M& m1, M1& m2) {
	M1 m3;
	m3.r = m1.r;
	m3.c = m2.c;
	//cout << m1.r << " " << m1.c <<" "<<m2.c<< endl;
	for (int i = 0; i < m1.r; i++) {
		for (int j = 0; j < m2.c; j++) {
			m3.ma[i][j] = 0;
			for (int k = 0; k < m1.c; k++) {
				m3.ma[i][j] = m3.ma[i][j] ^ (m1.ma[i][k] & m2.ma[k][j]);
			}
		}
	}
	//outputMatrix(m3);
	for (int i = 0; i < m2.r; i++) {
		for (int j = 0; j < m2.c; j++) {
			m2.ma[i][j] = m3.ma[i][j];
		}
	}
}
*/
void LowMC::clearMatrix(M& m) {
	for (int i = 0; i < m.r; i++) {
		for (int j = 0; j < m.c; j++) {
			m.ma[i][j] = 0;
		}
	}
}

void LowMC::clearMatrix1(M1& m) {
	for (int i = 0; i < m.r; i++) {
		for (int j = 0; j < m.c; j++) {
			m.ma[i][j] = 0;
		}
	}
}

void LowMC::matrixEq(M& src, M& dec) {
	dec.r = src.r;
	dec.c = src.c;
	for (int i = 0; i < src.r; i++) {
		for (int j = 0; j < src.c; j++) {
			dec.ma[i][j] = src.ma[i][j];
		}
	}
}

void LowMC::gauss(M& eqSys, int col) {
	int variableNum = col;
	bool isFirst = false;
	int targetRow = 0;

	for (int i = 0; i < variableNum; i++) {
		isFirst = true;
		for (int j = targetRow; j < eqSys.r; j++) {
			if (isFirst && eqSys.ma[j][i]) {
				isFirst = false;
				swap(eqSys.ma[j], eqSys.ma[targetRow]);
				targetRow++;
			}
			else {
				if (eqSys.ma[j][i]) {//apply Gauss
					for (int k = i; k < eqSys.c; k++) {
						eqSys.ma[j][k] ^= eqSys.ma[targetRow - 1][k];
					}
				}
			}
		}
	}
}

void LowMC::simplify(M& m, int col) {
	//further simplify
	int varNum = col;
	int index = 0;
	bool find = false;
	for (int i = 0; i < m.r; i++) {
		find = false;
		for (int j = index; j < varNum; j++) {
			if (m.ma[i][j]) {
				index = j;
				find = true;
				break;
			}
		}
		if (find) {//it is 1 in connectH[i][index], eliminate 1 above
			//cout << "find " << index << endl;
			for (int k = 0; k < i; k++) {
				if (m.ma[k][index]) {//add row[i] to row [k]
					for (int t = i; t < m.c; t++) {
						m.ma[k][t] = m.ma[k][t] ^ m.ma[i][t];
					}
				}
			}
			index++;
		}
	}
}

void LowMC::storeSolutions(vector<vector<bool> >& sol, M& eqSys, int& solNum) {
	vector<int> lead;
	vector<int> freebits;
	freebits.clear();
	lead.clear();
	bool* isFree;
	isFree = new bool[eqSys.c - 1];
	memset(isFree, 1, eqSys.c - 1);

	int start = 0;
	for (int r = 0; r < eqSys.r; r++) {
		while (start < eqSys.c - 1 && eqSys.ma[r][start] == 0) {
			start++;
		}
		if (start == eqSys.c - 1) {
			break;
		}
		lead.push_back(start);
		isFree[start] = false;
		start++;
	}

	if (lead.size() < eqSys.r) {
		for (int j = lead.size(); j < eqSys.r; j++) {
			if (eqSys.ma[j][eqSys.c - 1] != 0) {
				solNum = 0;
				return;
			}
		}
	}

	for (int i = 0; i < eqSys.c - 1; i++) {
		if (isFree[i]) {
			freebits.push_back(i);
		}
	}
	//cout << "free size:" << freebits.size() << endl;
	//cout << "lead size:" << lead.size() << endl;*/

	vector<bool> eachsol;
	eachsol.clear();
	eachsol.resize(eqSys.c - 1);
	int solSize = 1 << freebits.size();
	for (int i = 0; i < solSize; i++) {
		for (int j = 0; j < freebits.size(); j++) {
			eachsol[freebits[j]] = (i >> j) & 0x1;
		}
		for (int k = lead.size() - 1; k >= 0; k--) {
			//compute eachsol[lead[k]] use row= k
			eachsol[lead[k]] = eqSys.ma[k][eqSys.c - 1];
			for (int j = lead[k] + 1; j < eqSys.c - 1; j++) {
				if (eqSys.ma[k][j] == 1) {
					eachsol[lead[k]] = eachsol[lead[k]] ^ eachsol[j];
				}
			}
		}
		solNum++;
		sol.push_back(eachsol);
	}

	delete[]isFree;
	freebits.clear();
	lead.clear();
	eachsol.clear();
}

/*****************************
*
*output functions
*
*****************************/


void LowMC::outputM(M mat) {
	for (int i = 0; i < mat.r; i++) {
		for (int j = 0; j < mat.c; j++) {
			cout << mat.ma[i][j];
			//if (j == N - 1) {
				//cout << "  ";
			//}
		}
		cout << endl;
	}
}

void LowMC::outputArray(bool vec[], int size) {
	for (int i = 0; i < size; i++) {
		cout << vec[i];
	}
	cout << endl;
}

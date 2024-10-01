#ifndef _LOWMC_H_
#define _LOWMC_H_

#define N 129
#define P 43
#define R 4

const int D = N - 3 * P;

const int ddt[8][4] = {
	0,0,0,0,
	1,3,5,7,
	2,3,6,7,
	1,5,2,6,
	4,5,6,7,
	1,3,4,6,
	2,3,4,5,
	1,2,4,7
};

#include <string>
#include <vector>
using namespace std;


struct M {
	int r=N ;//default value
	int c=N ;//default value
	bool ma[400][400];
};
struct M1 {
	int r = N;//default value
	int c = N;//default value
	int ma[N][10000];
};
class LowMC {
private:
	vector<vector<bool> > CONS;
	M* L;
	M* IL;
	M* KF;
	M* IKF;
	//M WK;//whitening key
	M H1, Q1, GammaTr;
	M PP0, QQ1;
	vector<vector<unsigned long long>> DH;
public:
	LowMC();
	~LowMC();
	void initializePars(string filename, bool isKey, M* mat);

	//encryption
	void encryptionDetails(bool p0[], bool k[], bool c[], int rounds, vector<vector<bool> >& LOut, vector<vector<bool> >& SOut);

	
	//void constructForwardDiffR1(bool ga[]);
	void constructEquation2(vector<vector<bool> > SOutDiff, vector<vector<bool> > LOutDiff, int Eqnumber,int g1,int g2);

	//void constructEquation1(bool inputdiff[],int g1);

	void constructExpressions(bool diff[], M& eq,int col);

	void addCubic1(vector<vector<bool> >& qeq, int row, M& fac0, int in0, M& fac1, int in1, M& fac2, int in2, int total,int col);

	void addLinear1(vector<vector<bool> >& qeq, int row, M& lin, int index,int col);

	void addQuadratic1(vector<vector<bool> >& qeq, int row, M& fac0, int in0, M& fac1, int in1, int total,int col);

	//void addCubic2(bool  qeq[][1351], int row, M& fac0, int in0, M& fac1, int in1, M& fac2, int in2, int total);

	//void addLinear2(vector<vector<bool> >& qeq, int row, M& lin, int index, int col);
	//void addQuadratic2(bool  qeq[][1351], int row, M& fac0, int in0, M& fac1, int in1, int total);
	//bool checkCorrectness(vector<bool>& gamma, vector<bool>& alpha, vector<bool>& u,int rounds,int total);

	//key recovery
	void keyRecovery(bool key[], bool c[], int r, vector<vector<bool> >& LOut, vector<vector<bool> >& SOut, vector<vector<bool> >& sOut, vector<vector<bool> >& lOut, bool y[]);
	void dynamicallyUpdateExpression(M& expLOut, M& expSOut, int in, int out,int Sbox);
	void extractEquations(M& eq, M& expSOut, int in, int out,int Sbox);
	void solveKey(M& eqSys, M& extraEq, int extraSize, bool testVec[], bool key[], bool yy[]);
	void addLinear(M& qeq, int row, M& lin, int index);
	void addQuadratic(M& qeq, int row, M& fac0, int in0, M& fac1, int in1, int total);

	//search algebraic equations
	void searchAlgebraicEquations();
	void computeQuadraticEq(bool x[], bool row[]);
	void computeCubicEq(bool x[], bool row[]);

	//matrix operations
	void matrixMul(M& m, bool x[], bool y[]);
	void matrixMul(M& m, vector<bool>& x, vector<bool>& y);
	void matrixMul(M& m, bool x[]);
	void matrixMul(M& m1, M& m2, M& m3);
	void matrixMul(M& m1, M& m2);
	void matrixMul2(M& m1, M1& m2);
	void clearMatrix(M& m);
	void clearMatrix1(M1& m);
	void matrixEq(M& src, M& dec);
	void gauss(M& eqSys, int col);
	void simplify(M& eqSys, int col);
	void storeSolutions(vector<vector<bool> >& sol, M& eqSys, int& solNum);

	//output functions
	void outputM(M mat);
	void outputArray(bool vec[], int size);
};

#endif

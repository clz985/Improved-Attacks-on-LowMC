#include "LowMC.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <ctime>
#include <random>
using namespace std;



void testSuccessRate() {
	LowMC lowmc;


	bool p0[N], k[N], c0[N], p1[N], c1[N];
	vector<vector<bool> > SOut0, SOut1, LOut0, LOut1;
	SOut0.resize(R), SOut1.resize(R), LOut0.resize(R), LOut1.resize(R);
	for (int i = 0; i < R; i++) {
		SOut0[i].resize(N), SOut1[i].resize(N), LOut0[i].resize(N), LOut1[i].resize(N);
	}
	vector<vector<bool> > SOutDiff, LOutDiff;
	SOutDiff.resize(R), LOutDiff.resize(R);
	for (int i = 0; i < R; i++) {
		SOutDiff[i].resize(N), LOutDiff[i].resize(N);
	}

	srand(time(NULL));
	int success = 0;
	int CNT = 800;

	std::mt19937 mt;
	std::random_device rnd;
	mt.seed(rnd());

	for (int cnt = 0; cnt < CNT; cnt++) {
		for (int i = 0; i < N / 3; i++) {
			unsigned int v0 = 0, v1 = 0;
			v0 = mt(), v1 = mt();
			//v0 = 0xffff, v1 = 0xffff;
			//v0 = 0xffff, v1 = 0;
			for (int j = 0; j < 3; j++) {
				p0[i * 3 + j] = (v0 >> j) & 0x1;
				k[i * 3 + j] = (v1 >> j) & 0x1;
			}
		}

		lowmc.encryptionDetails(p0, k, c0, R, LOut0, SOut0);

		bool diff[N] = {
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
		};
		
		for (int i = 0; i < N; i++) {
			p1[i] = p0[i] ^ diff[i];
		}
		lowmc.encryptionDetails(p1, k, c1, R, LOut1, SOut1);



		int rounds = R;
		int inactiveCnt[3] = { 0,0,0 };

		for (int i = rounds; i > 1; i--) {
			for (int j = 0; j < N / 3; j++)
			{
				if (SOut0[i - 1][3 * j] == SOut1[i - 1][3 * j]
					&& SOut0[i - 1][3 * j + 1] == SOut1[i - 1][3 * j + 1]
					&& SOut0[i - 1][3 * j + 2] == SOut1[i - 1][3 * j + 2]) {
					inactiveCnt[rounds - i]++;
				}
			}
	
		}

		if (inactiveCnt[0] <= 11 && inactiveCnt[1] <= 11 && inactiveCnt[2] <= 11) {
			if (inactiveCnt[0] >= 4 && inactiveCnt[2] >= 4)
			{

				success++;

				for (int i = 0; i < R; i++) {
					for (int j = 0; j < N; j++) {
						SOutDiff[i][j] = SOut0[i][j] ^ SOut1[i][j];
						LOutDiff[i][j] = LOut0[i][j] ^ LOut1[i][j];
					}
				}

			}
		}



		/*if (i % 100 == 0) {
			cout << i << ":" << success << endl;
		}*/
	}

	cout << "rate:" << (double)(success) / CNT << endl;

	for (int i = 0; i < R; i++) {
		SOut0[i].clear(), SOut1[i].clear(), LOut0[i].clear(), LOut1[i].clear();
	}
	SOut0.clear(), SOut1.clear(), LOut0.clear(), LOut1.clear();
}


void constructEq2() {
	LowMC lowmc;
	

	bool p0[N], k[N], c0[N], p1[N], c1[N];
	vector<vector<bool> > SOut0, SOut1, LOut0, LOut1;
	SOut0.resize(R), SOut1.resize(R), LOut0.resize(R), LOut1.resize(R);
	for (int i = 0; i < R; i++) {
		SOut0[i].resize(N), SOut1[i].resize(N), LOut0[i].resize(N), LOut1[i].resize(N);
	}
	vector<vector<bool> > SOutDiff, LOutDiff;
	SOutDiff.resize(R), LOutDiff.resize(R);
	for (int i = 0; i < R; i++) {
		SOutDiff[i].resize(N), LOutDiff[i].resize(N);
	}

	srand(time(NULL));
	int success = 0;


	std::mt19937 mt;
	std::random_device rnd;
	mt.seed(rnd());

	
		for (int i = 0; i < N / 3; i++) {
			unsigned int v0 = 0, v1 = 0;
			v0 = mt(), v1 = mt();
			//v0 = 0xffff, v1 = 0xffff;
			//v0 = 0xffff, v1 = 0;
			for (int j = 0; j < 3; j++) {
				p0[i * 3 + j] = (v0 >> j) & 0x1;
				k[i * 3 + j] = (v1 >> j) & 0x1;
			}
		}

		lowmc.encryptionDetails(p0, k, c0, R, LOut0, SOut0);
		
		bool diff[N] = {
			0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,
		};
	
		for (int i = 0; i < N; i++) {
			p1[i] = p0[i] ^ diff[i];
		}
		lowmc.encryptionDetails(p1, k, c1, R, LOut1, SOut1);

	
		int rounds = R;
		
		//int inactiveCnt = 0;
		//check the last 98 S-boxes
		
			for (int i = 0; i < R; i++) {
				for (int j = 0; j < N; j++) {
					SOutDiff[i][j] = SOut0[i][j] ^ SOut1[i][j];
					LOutDiff[i][j] = LOut0[i][j] ^ LOut1[i][j];
				}
			}
			/*
			cout << " After S-box" << endl;
			for (int i = 0; i < R; i++) {
				for (int j = 0; j < N; j++) {
					cout<<SOutDiff[i][j];
					
				}
				cout << endl;
			}
			cout << " After Linear" << endl;
			for (int i = 0; i < R; i++) {
				for (int j = 0; j < N; j++) {
					cout << LOutDiff[i][j];

				}
				cout << endl;
			}
			*/
			lowmc.constructEquation2(SOutDiff, LOutDiff, 1, 35, 13);

	cout << "find" << endl;

	for (int i = 0; i < R; i++) {
		SOut0[i].clear(), SOut1[i].clear(), LOut0[i].clear(), LOut1[i].clear();
	}
	SOut0.clear(), SOut1.clear(), LOut0.clear(), LOut1.clear();
}
int main() {
	//cout << "1 -> test the construction with a full S-box layer (800 tests in a few minutes)" << endl;
	//cout << "2 -> construct the cubic equation system" << endl;
	int cmd;
	cout << endl << "input command (1/2):";
	cin >> cmd;
	if(cmd==1)
	{
	testSuccessRate();
	}
	if (cmd == 2)
	{
		constructEq2();
	}
	return 0;
}
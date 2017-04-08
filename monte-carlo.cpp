#include<iostream>
#include<time.h>
#include<cstdarg>
#include<random>
#include<map>
#include<fstream>
#include<string>
#include <math.h>
#include<iomanip>

using namespace std;
#define N 5
#define NDIM 2
#define NS pow(N,NDIM)
#define N2 NS/2
#define NCONFIG 2
#define NSQ2 NS/2
#define NTEMP 2


float uniformdist(void) {
	random_device rd;
	mt19937 gen(rd());

	// values near the mean are the most likely
	// standard deviation affects the dispersion of generated values from the mean
	uniform_real_distribution<> d(0,1);
	
	return d(gen);
	/*map<int, int> hist;
	for (int n = 0; n<10000; ++n) {
		++hist[(d(gen))];
	}
	for (auto p : hist) {
		std::cout << std::fixed << std::setprecision(1) << std::setw(2)
			<< p.first << ' ' << '\n';
	}*/
}
long double myrand(long iseed) {
	return uniformdist();
	

}

void write(int fname, int *i) {
	ofstream f;
	f.open(to_string(fname)+".txt", ios::out);
	for (int j = 0; j < 5; j++) {
		f << *(i+j) <<endl ;
		cout << "writing  " << *(i+j) << endl;
		
	}
	
	f.close();
}
void printmat(float (*phi)[N ]) {
	for (int i = 0; i < N ; i++) {
		for (int j = 0; j < N ; j++) {
			cout << phi[i][j] << "\t";
		}
		cout << endl;
	}
}
void printmat1(int(*occ)[N ]) {
	for (int i = 0; i < N ; i++) {
		for (int j = 0; j < N ; j++) {
			cout << occ[i][j] << "\t";
		}
		cout << endl;
	}
}

///  Initializing phi[][]  randomly
void rsiteen(int n, int iseed,float (*phi)[N] ,float w) {
	
	
	int i, j, n2, k;
	float aa;
	for (i = 1; i < N-1; i++) {
		for (j = 1; j < N-1; j++) {
			aa = myrand( iseed);
			//cout << aa << endl;
			if (aa <= 0.5) {
				phi[i][j] = -1 * w;
			}
			else {
				phi[i][j] = w;
			}
			
		}
	}
	printmat(phi);
	

}
///// Initializing occ[][] randomly
void rocc(int n,int iseed,int (*occ)[N]) {
	float a;
	
	
	for (int i = 1; i < N-1; i++){
		for (int j = 1; j < N-1; j++) {
			a = 2 * myrand(iseed)-1;
			if (a > 0) {
				occ[i][j] = -1;
			}
			else {
				occ[i][j] = 1;
			}
		}
	}

	for (int j = 1; j < N-1; j++) {
		occ[N-1][j] = occ[1][j];
		occ[0][j] = occ[N-2][j];
	}
	for (int i = 1; i < N-1; i++) {
		occ[i][N-1] = occ[i][1];
		occ[i][0] = occ[i][N-2];
	}

	
}

void hamiltonian(int n, int occ[N][N],float phi[N][N],float &toten ,float &totenst) {
	float toten1;
	toten = 0;
	for (int i = 1; i < n-1 ; i++) {
		for (int j = 1; j < n-1 ; j++) {
			toten = toten - phi[i][j] * occ[i][j];
		}
	}
	toten1 = 0;
	for (int ix =1 ; ix < n-1 ; ix++) {
		for (int iy = 0; iy < n-1 ; iy++) {
			toten1 = toten1 - occ[ix][iy] * (occ[(ix) + 1][iy] + occ[(ix) - 1][iy] + occ[ix][(iy) + 1] + occ[ix][(iy) - 1]);
		}
	}
	toten = toten + toten1 / 2;
	totenst = toten / n*n;
	cout << "totenst" << totenst << endl;

}


///  Stores site energy which corresponds to 4 neighbours in    sten1[][]
void siteEnergy(int n, int occ[N][N], float phi[N][N], float (*sten1)[N]) {
	for (int ix = 1; ix < n - 1; ix++) {
		for (int iy = 1; iy < n - 1; iy++) {
			sten1[ix][iy] = -1 * phi[ix][iy];
			sten1[ix][iy] = sten1[ix][iy] - occ[(ix) + 1][iy] - occ[(ix) - 1][iy] - occ[ix][(iy) - 1] - occ[ix][(iy) + 1];

		}
	}
	printmat(sten1);
	cout << endl;
}

void updatehe(int n, int ic, int jc, float (*sten1)[N], int occ[N ][N]) {
	//sten1[i]
	sten1[ic - 1][jc] += 2 * occ[ic][jc];
	sten1[ic + 1][jc] += 2 * occ[ic][jc];
	sten1[ic][jc - 1] += 2 * occ[ic][jc];
	sten1[ic][jc + 1] += 2 * occ[ic][jc];
	//printmat(sten1);
	//cout << endl;

}

int main1(void) {


	int iseed, x1, ix, ic, iy, jc, imeas, iskip, nskip, nmcs, k, nmeas;
	int it, itemp, ntime, occ[N ][N ] = { 0 }, kconfig, nmeass[NTEMP];
	float phi[N][N] = { 0 }, aaa, tempp[NTEMP] = { 0 }, phisum;
	float toten, totenst, beta, a, c, delen, ran2, mag, temp;
	float sten1[N][N] = { 0 }, delst, start, finish, totenh;
	ofstream f2,f3,f4,f200,f300,f400,f12;
	f2.open("mylog/2.txt", ofstream::app);
	f2 << "i" << "\t" << "j" << "\t" << "phi[i][j]" << endl;
	f3.open("mylog/3.txt", ofstream::app);
	f3 << "i" << "\t" << "j" << "\t" << "occ[i][j]" << endl;
	f4.open("mylog/4.txt", ofstream::app);
	int flip = 0;
	float w = 0.0;
	iseed = 5;
	
	nskip = 1000;
	for (itemp = 0; itemp < NTEMP; itemp++) {
		tempp[itemp] = 3.5 - itemp;
		nmeass[itemp] = 100;
		

	}

	for (kconfig = 0; kconfig < NCONFIG; kconfig++) {
		cout << "config " << kconfig << endl;
		rsiteen(N, iseed, phi, w); 
		cout << w << endl;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				
				f2 << i <<"\t" << j << "\t" <<phi[i][j] << endl;
			}
		}
		f2 << endl;
		phisum = 0;
		for (int i = 0; i < N ; i++) {
			for (int j = 0; j < N ; j++) {
				phisum += phi[i][j];
			}
		}

		rocc(N, iseed, occ);
		for (int i = 0; i < N ; i++) {
			for (int j = 0; j < N ; j++) {
				f3 << i << "\t" << j << "\t" << occ[i][j] << endl;

			}
		}
		f3 << endl;
		siteEnergy(N, occ, phi, sten1);
		hamiltonian(N, occ, phi, toten, totenst);
		cout << totenst << "main" << endl;
		f4 << "initialen-toten" << "\t" << toten << "\t" <<"encrystal=" << "\t" << 2 * N*N << "\t" << phisum << "\t"  <<"endisordercrystal=" <<"\t" << - 2 * N*N-phisum  <<"\t"<<   -2 * N*N + phisum << endl;
		f4 << endl;
		printmat1(occ);
		cout << endl;
		for (itemp = 0; itemp < NTEMP; itemp++) {
			cout << "temp : " << tempp[itemp] << endl;
			cout << "No of Simul : " << nmeass[itemp] << endl;
			nmeas = nmeass[itemp];
			temp = tempp[itemp];
			beta = 1 / temp;
			toten = 0;
			/***************************        Monte Carlo Simulation      *****************************/
			for (imeas = 0; imeas < nmeas; imeas++) {
				//cout <<  "  itemp : " << itemp << "  imeas :" << imeas << endl;
				cout << "imeas :" << imeas << endl;
				for (iskip = 0; iskip < nskip; iskip++) {
					//printmat1(occ);
					//cout << "iskip : " << iskip << endl;
					//cout << "temp : " << nmeass[itemp] << "  itemp : " << itemp  << "  imeas :" << imeas << "  iskip : " << iskip <<  endl;
					for (nmcs = 0; nmcs < NS; nmcs++) {
						//printmat1(occ);
						ic = int(myrand(iseed)*(N-1));
						jc = int(myrand(iseed)*(N-1));
						if ((ic == 0) || (ic == N-1) || (jc == 0) || (jc == N-1)) {
							continue;
						}
						delst = -2 * occ[ic][jc] * sten1[ic][jc];

						if (delst <= 0) {
							toten += delst;
							occ[ic][jc] = -1*occ[ic][jc];
							updatehe(N, ic, jc, sten1, occ);
							flip++;

						}
						else {
							a = exp(-delst*beta);
							c = myrand(iseed);
							if (c <= a) {
								toten += delst;
								occ[ic][jc] = -1*occ[ic][jc];
								updatehe(N, ic, jc, sten1, occ); 
								flip++;
								
							}
						}
					}
				}
				mag = 0;
				for (int i = 0; i < N ; i++) {
					for (int j = 0; j < N ; j++) {
						mag += occ[i][j];
					}
				}

				f400.open("mylog/40" + to_string(kconfig) + ".txt", ofstream::out);
				f400 << "temp" << "\t" << "imeas" << "\t" << "mag / NS" << "\t" << "toten" << endl;
				f400 << temp << "\t" << imeas << "\t" << mag / NS << "\t" << toten << endl;
				f400 << endl;
				f400.close();
			}

			mag = 0;
			for (int i = 0; i < N ; i++) {
				for (int j = 0; j < N; j++) {
					mag += occ[i][j];
				}
			}
			f4 << "kconfig" << "\t" << "temp" << "\t" << "\t" << "mag / NS" << "\t" << "toten" << endl;
			f4 << kconfig << "\t" << temp << "\t" << "\t" << mag / NS << "\t" << toten << endl;
		}
		//int clusters = hoshen_kopelman(matrix, N,N);
		
		f200.open("mylog/20" + to_string(kconfig) + ".txt", ofstream::out);
		f300.open("mylog/30" + to_string(kconfig) + ".txt", ofstream::out);
		for (int i = 0; i < N ; i++) {
			for (int j = 0; j < N ; j++) {
				if (occ[i][j] < 0) {
					
					f200 << i << "\t" << j << "\t" << occ[i][j] << endl;
				}
				else {
					
					f300 << i << "\t" << j << "\t" << occ[i][j] << endl;
				}
			}
		}
		f200.close();
		f300.close();
	
		f12 << kconfig << mag / NS << toten << endl;
		printmat1(occ);

	}
	f12.close();
	f2.close();
	f3.close();
	//printmat(phi);
	cout << endl;
	printmat(sten1);
	//printmat1(occ);
	return 0;

}
void printz(float(*z)[2]) {
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			//z[i][j] += 10;
			cout << z[i][j] << endl;
		} 
	}
}
void main(void) {
	
	main1();
	
	//a(ab);
	//printz(ab);*/
	system("pause");
}

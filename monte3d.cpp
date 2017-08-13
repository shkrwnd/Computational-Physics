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
#define tab			"\t"

#define N 10
#define NDIM 3
#define NS pow(N,NDIM)
#define N2 NS/2
#define NCONFIG 2
#define NSQ2 NS/2
#define NTEMP 6
#define NMEAS 2+1
#define INFINITE 10000
#define HIGH INFINITE
#define sqlat0 occ
#define sqlat1 occ
#define VER N-1
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))



float uniformdist(void) {
	random_device rd;
	mt19937 gen(rd());
	// values near the mean are the most likely
	// standard deviation affects the dispersion of generated values from the mean
	uniform_real_distribution<> d(0,1);
	return d(gen);
	
}
long double myrand(long iseed) {
	return uniformdist();
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
void rsiteen(int n, int iseed,float (*phi)[N][N] ,float w) {
	
	
	int i, j, n2, k;
	float aa;
	for (i = 1; i < N-1; i++) {
		for (j = 1; j < N-1; j++) {
			for (int k = 1;k< N-1;k++){
				aa = myrand( iseed);

				if (aa <= 0.5) {
					phi[i][j][k] = -1 * w;
				}
				else {
					phi[i][j][k] = w;
				}
			}
			
		}
	}
	//printmat(phi);
	

}

///// Initializing occ[][] randomly
void rocc(int n,int iseed,int (*occ)[N][N]) {
	float a;
	
	
	for (int i = 1; i < N-1; i++){
		for (int j = 1; j < N-1; j++) {
			for (int k = 1;k < N-1 ;k++){
				a = 2 * myrand(iseed)-1;
				if (a > 0) {
					occ[i][j][k] = -1;
				}
				else {
					occ[i][j][k] = +1;
				}
			}
		}
	}

	for (int k = 1 ;k < N-1 ;k++){
		for (int j = 1; j < N-1; j++) {
			occ[N-1][j][k] = occ[1][j][k];
			occ[0][j][k] = occ[N-2][j][k];

		}
		for (int i = 1; i < N-1; i++) {
			occ[i][N-1][k] = occ[i][1][k];
			occ[i][0] = occ[i][N-2][k];
		}
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
		for (int iy = 1; iy < n-1 ; iy++) {
			toten1 = toten1 - occ[ix][iy] * (occ[(ix) + 1][iy] + occ[(ix) - 1][iy] + occ[ix][(iy) + 1] + occ[ix][(iy) - 1]);
		}
	}
	toten = toten + toten1 / 2;
	totenst = toten / n*n;
	//cout << "toten" << toten << endl;

}


///  Stores site energy which corresponds to 4 neighbours in    sten1[][]
void siteEnergy(int n, int occ[N][N], float phi[N][N], float (*sten1)[N]) {
	for (int ix = 1; ix < n - 1; ix++) {
		for (int iy = 1; iy < n - 1; iy++) {
			for(int iz = 1;iz < n-;iz++){
				sten1[ix][iy][iz] = -1 * phi[ix][iy];
				sten1[ix][iy] = sten1[ix][iy] - (occ[(ix) + 1][iy] + occ[(ix) - 1][iy] + occ[ix][(iy) - 1] + occ[ix][(iy) + 1]);
			}
			/*cout << sten1[ix][iy] << endl;
			cout << phi[ix][iy] << endl;*/
		}
	}
	/*cout << "occ matt" << endl;
	printmat1(occ);
	cout << "sten mat" << endl;
	printmat(sten1);
	cout << endl;*/
	
	
}

void updatehe(int n, int ic, int jc, float (*sten1)[N], int occ[N ][N]) {
	//sten1[i]
	if(ic-1 !=0)
		sten1[ic - 1][jc] += 2 * occ[ic][jc];
	if(ic+1 !=N-1)
		sten1[ic + 1][jc] += 2 * occ[ic][jc];
	if(jc-1 != 0)
		sten1[ic][jc - 1] += 2 * occ[ic][jc];
	if(jc+1 != N-1)
		sten1[ic][jc + 1] += 2 * occ[ic][jc];
	//printmat(sten1);
	//cout << endl;

}

int main1(void) {


	int iseed, x1, ix, ic, iy, jc, imeas, iskip, nskip, nmcs, k, nmeas;
	int it, itemp, ntime, occ[N ][N ][N] = { 0 }, kconfig, nmeass[NTEMP];
	float phi[N][N][N] = { 0 }, aaa, tempp[NTEMP] = { 0 }, phisum;
	float toten, totenst, beta, a, c, delen, ran2, mag, temp;
	float sten1[N][N][N] = { 0 }, delst, start, finish, totenh;
	ofstream f2,f3,f4,f200,f300,f400,f12;
	f2.open("mylog/2.txt", ofstream::app);
	f2 << "i" << "\t" << "j" << "\t" << "phi[i][j]" << endl;
	f3.open("mylog/3.txt", ofstream::app);
	f3 << "i" << "\t" << "j" << "\t" << "occ[i][j]" << endl;
	f4.open("mylog/4.txt", ofstream::app);
	int flip = 0;
	float w = 0.0;
	iseed = 5;
	int mag2 = 0;
	nskip = 1000;
	for (itemp = 0; itemp < NTEMP; itemp++) {
		tempp[itemp] =3.5 - itemp*0.5;
		nmeass[itemp] = NMEAS;
		 

	}

	for (kconfig = 0; kconfig < NCONFIG; kconfig++) {
		cout << "config " << kconfig << endl;
		rsiteen(N, iseed, phi, w); 
		cout << w << endl;
		/*for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				
				f2 << i <<"\t" << j << "\t" <<phi[i][j] << endl;
			}
		}*/
		f2 << endl;
		phisum = 0;
		/*for (int i = 0; i < N ; i++) {
			for (int j = 0; j < N ; j++) {
				phisum += phi[i][j];
			}
		}*/

		rocc(N, iseed, occ);
		/*for (int i = 0; i < N ; i++) {
			for (int j = 0; j < N ; j++) {
				f3 << i << "\t" << j << "\t" << occ[i][j] << endl;

			}
		}*/
		f3 << endl;
		siteEnergy(N, occ, phi, sten1);
		hamiltonian(N, occ, phi, toten, totenst);
		//cout << toten << "main" << endl;
		f4 << "initialen-toten" << "\t" << toten << "\t" <<"encrystal=" << "\t" << 2 * N*N << "\t" << phisum << "\t"  <<"endisordercrystal=" <<"\t" << - 2 * N*N-phisum  <<"\t"<<   -2 * N*N + phisum << endl;
		f4 << endl;
		//printmat1(occ);
		cout << endl;
		for (itemp = 0; itemp < NTEMP; itemp++) {
			cout << "temp : " << tempp[itemp] << endl;
			cout << "No of Simul : " << nmeass[itemp] << endl;
			nmeas = nmeass[itemp];
			temp = tempp[itemp];
			beta = 1 / temp;
			//toten = 0;
			/***************************        Monte Carlo Simulation      *****************************/
			for (imeas = 0; imeas < nmeas; imeas++) {

				//cout <<  "  itemp : " << itemp << "  imeas :" << imeas << endl;
				cout << "imeas :" << imeas << endl;
				//printmat1(occ);
				cout << endl;
				//printmat(sten1);
				cout  << "toten" << "\t" << toten  << endl;
				//cout << "phisum" << "\t" << phisum << endl;
				for (iskip = 0; iskip < nskip; iskip++) {
					//printmat1(occ);
					//cout << "iskip : " << iskip << endl;
					//cout << "temp : " << nmeass[itemp] << "  itemp : " << itemp  << "  imeas :" << imeas << "  iskip : " << iskip <<  endl;
					for (nmcs = 0; nmcs < NS; nmcs++) {
						//printmat1(occ);
						ic = int(myrand(iseed)*(N-1));
						jc = int(myrand(iseed)*(N-1));
						kc = int(myrand(iseed)*(N-1));
						if ((ic != 0) && (ic != N-1) && (jc != 0) && (jc != N-1) && (kc != N-1) && (kc != 0)) {
							//delst = -2 * occ[ic][jc] * sten1[ic][jc];
							//delst = -2 * occ[ic][jc] ;
							delst =  2 * occ[ic][jc] * (occ[ic + 1][jc] + occ[ic - 1][jc] + occ[ic][jc + 1] + occ[ic][jc - 1]) - 2*phi[ic][jc]*occ[ic][jc];
							if (delst <= 0) {
								toten += delst;
								occ[ic][jc] = -1 * occ[ic][jc];
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
				}
				mag = 0;
				for (int i = 0; i < N ; i++) {
					for (int j = 0; j < N ; j++) {
						mag += occ[i][j];
					}
				}
				
				cout << mag << "\t" <<  "mag" << endl;
				cout << mag2 << "\t" << "del(m)" << endl;
				mag2 = mag;
				
				f400.open("mylog/40" + to_string(kconfig) + ".txt", ofstream::app);
				f400 << "temp" << "\t" << "imeas" << "\t" << "mag / NS" << "\t" << "toten" << endl;
				f400 << temp << "\t" << imeas << "\t" << mag / NS << "\t" << toten << endl;
				f400 << endl;
				f400.close();
			}
			//cluster call
			//clustering(occ);
			mag = 0;
			for (int i = 0; i < N ; i++) {
				for (int j = 0; j < N; j++) {
					mag += occ[i][j];
				}
			}
			//cout << mag << endl;
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
		

	}
	f12.close();
	f2.close();
	f3.close();
	//printmat(phi);
	cout << flip << "\t" << "flips" << endl;
	//printmat(sten1);
	//printmat1(occ);
	return 0;

}
void printz(float(*z)[N]) {
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			//z[i][j] += 10;
			cout << z[i][j] << endl;
		} 
	}
}
void main(void) {
	
	int a[N][N][N] = {0};
	rocc(N,5,a);
	
	//main1();
	system("pause");
}

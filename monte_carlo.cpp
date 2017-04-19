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
#define NDIM 2
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

int gotoindex(int A[N - 1], int index)
{
	int k = 0;
	if (A[index] < 0)
	{
		//cout << "\nconflict of" << A[index] << "index changed to" << abs(A[index]) << endl;
		gotoindex(A, abs(A[index]));
	}
	else
	{
		//cout << "\nindex ret: " << index << endl;
		return index;
	}
}
void clustering(int occ[N][N]) {
	int clusstats1[N / 2] = { 0 };
	int clusstats0[N / 2] = { 0 };
	//agumented matrix
	int i = 0, j = 0;
	for (i = 0; i <= VER; i++)
	{

		sqlat0[i][0] = HIGH;
		sqlat1[i][0] = HIGH;
	}
	for (j = 0; j <= VER; j++)
	{

		sqlat0[0][j] = HIGH;
		sqlat1[0][j] = HIGH;
	}

	for (i = 1; i <= VER; i++)
	{
		for (j = 1; j <= VER; j++)
		{
			//sqlat0[i][j] = visited[(i - 1)*VER + (j - 1) + 1];
			//sqlat1[i][j] = visited[(i - 1)*VER + (j - 1) + 1];

			if (sqlat1[i][j] == 0)sqlat1[i][j] = HIGH;
			if (sqlat0[i][j] == 1)sqlat0[i][j] = HIGH;

			//cout << visited[(i - 1)*VER + (j - 1) + 1];

		}
	}
	//cout << "\nsqlat 0\n";
	//for (i = 0; i <= VER; i++)
	//{
	//	for (j = 0; j <= VER; j++)
	//	{
	//		if (sqlat0[i][j] == HIGH)
	//			//cout << "- ";
	//			cout << sqlat0[i][j] << " ";
	//		else
	//			cout << sqlat0[i][j] << " ";
	//			//cout << "  ";
	//	}cout << endl;
	//}
	//cout << "\nsqlat 1\n";
	//for (i = 0; i <= VER; i++)
	//{
	//	for (j = 0; j <= VER; j++)
	//	{
	//		if (sqlat1[i][j] == HIGH)
	//			//cout << "- ";
	//			cout << sqlat1[i][j] << " ";
	//		else
	//			cout << sqlat1[i][j] << " ";
	//		//cout << "  ";
	//	}cout << endl;
	//}
	//clustering
	//cluster 1
	int m = 0, n = 0, min = 0, a = 0, b = 0;

	for (i = 1; i < VER + 1; i++)
	{

		for (j = 1; j < VER + 1; j++)
		{

			if (sqlat1[i][j] == 1)
			{
				a = sqlat1[i - 1][j];
				b = sqlat1[i][j - 1];
				//cout << "\n\nbw " << a << "(" << i - 1 << "," << j << ") & " << b << " (" << i << "," << j - 1 << " )\n";
				min = MIN(a, b);
				if (a == HIGH && b == HIGH)
				{//rule 1
				 //cout << "rule:" << 1 << d;
					m++;
					clusstats1[m]++;
					sqlat1[i][j] = m;
					//cout << "ret:" << m << d;
				}
				else if ((a == HIGH && b < HIGH) || (a < HIGH && b == HIGH))
				{//rule 2
				 //cout << "rule:" << 2 << d;
					min = MIN(a, b);
					n = gotoindex(clusstats1, min);
					clusstats1[n]++;
					sqlat1[i][j] = n;
				}
				else if (a < b)
				{//rule 3
				 //cout << "rule:" << 3 << d;
					n = gotoindex(clusstats1, a);
					if (a != gotoindex(clusstats1, b))
						clusstats1[n] += clusstats0[gotoindex(clusstats1, b)] + 1;
					else
						clusstats1[n] += 1;
					clusstats1[b] = -n;
					sqlat1[i][j] = n;
				}
				else if (a == b)
				{//rule 4
				 //cout << "rule:" << 4 << d;
					n = gotoindex(clusstats1, a);
					clusstats1[n] += 1;
					sqlat1[i][j] = n;
				}
				else
				{//rule 5
				 //cout << "rule:" << 5 << d;
					n = gotoindex(clusstats1, b);
					if (b != gotoindex(clusstats1, a))
						clusstats1[n] += clusstats0[gotoindex(clusstats1, a)] + 1;
					else
						clusstats1[n] += 1;
					clusstats1[a] = -n;
					sqlat1[i][j] = n;
				}
				/*cout << "\ncluster\n\n";

				for (int ki = 0; ki <= VER; ki++)
				{
				for (int kj = 0; kj <= VER; kj++)
				{
				if (sqlat[ki][kj] == HIGH)
				cout << "- ";

				else
				cout << sqlat[ki][kj] << " ";
				}cout << endl;
				}*/
				/*cout << "clus stats\n";
				for (int k = 0; k <= m; k++)
				{
				cout << clusstats[k] << tab;
				}*/
			}
		}

	}

	//cluster 0
	m = 0, n = 0, min = 0, a = 0, b = 0;
	for (i = 1; i < VER + 1; i++)
	{

		for (j = 1; j < VER + 1; j++)
		{

			if (sqlat0[i][j] == 0)
			{
				a = sqlat0[i - 1][j];
				b = sqlat0[i][j - 1];
				//cout << "\n\nbw " << a << "(" << i - 1 << "," << j << ") & " << b << " (" << i << "," << j - 1 << " )\n";
				min = MIN(a, b);
				if (a == HIGH && b == HIGH)
				{//rule 1
				 //cout << "rule:" << 1 << d;
					m++;
					clusstats0[m]++;
					sqlat0[i][j] = m;
					//cout << "ret:" << m << d;
				}
				else if ((a == HIGH && b < HIGH) || (a < HIGH && b == HIGH))
				{//rule 2
				 //cout << "rule:" << 2 << d;
					min = MIN(a, b);
					n = gotoindex(clusstats0, min);
					clusstats0[n]++;
					sqlat0[i][j] = n;
				}
				else if (a < b)
				{//rule 3
				 //cout << "rule:" << 3 << d;
					n = gotoindex(clusstats0, a);
					if (a != gotoindex(clusstats0, b))
						clusstats0[n] += clusstats0[gotoindex(clusstats0, b)] + 1;
					else clusstats0[n] += 1;
					clusstats0[b] = -n;
					sqlat0[i][j] = n;
				}
				else if (a == b)
				{//rule 4
				 //cout << "rule:" << 4 << d;
					n = gotoindex(clusstats0, a);
					clusstats0[n] += 1;
					sqlat0[i][j] = n;
				}
				else
				{//rule 5
				 //cout << "rule:" << 5 << d;
					n = gotoindex(clusstats0, b);
					if (b != gotoindex(clusstats0, a))
						clusstats0[n] += clusstats0[gotoindex(clusstats0, a)] + 1;
					else clusstats0[n] += 1;
					clusstats0[a] = -n;
					sqlat0[i][j] = n;
				}
				/*cout << "\ncluster\n\n";

				for (int ki = 0; ki <= VER; ki++)
				{
				for (int kj = 0; kj <= VER; kj++)
				{
				if (sqlat[ki][kj] == HIGH)
				cout << "- ";

				else
				cout << sqlat[ki][kj] << " ";
				}cout << endl;
				}*/
				/*cout << "clus stats\n";
				for (int k = 0; k <= m; k++)
				{
				cout << clusstats[k] << tab;
				}*/
			}
		}

	}
	cout << "\ncluster\n\n";

	//for (i = 0; i <= VER; i++)
	//{
	//	for (j = 0; j <= VER; j++)
	//	{
	//		if (sqlat[i][j] == HIGH)
	//			cout << "- ";

	//		else
	//			cout << sqlat[i][j] << " ";
	//			//cout <<"  ";
	//	}cout << endl;
	//}
	//
	cout << "clus stats for 0\n";

	for (i = 0; i < N; i++)
	{

		cout << clusstats0[i] << tab;
	}cout << "\n\n";
	cout << "clus stats for 1\n";

	for (i = 0; i < N; i++)
	{

		cout << clusstats1[i] << tab;
	}cout << "\n\n";

}

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
				occ[i][j] = +1;
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
			sten1[ix][iy] = -1 * phi[ix][iy];
			sten1[ix][iy] = sten1[ix][iy] - (occ[(ix) + 1][iy] + occ[(ix) - 1][iy] + occ[ix][(iy) - 1] + occ[ix][(iy) + 1]);
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
				printmat1(occ);
				cout << endl;
				printmat(sten1);
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
						if ((ic != 0) && (ic != N-1) && (jc != 0) && (jc != N-1)) {
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
	
	main1();
	system("pause");
}

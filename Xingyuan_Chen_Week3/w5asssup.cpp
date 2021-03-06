#include <iostream>
#include <mpi.h>
#include <math.h>  
#include <vector>
using namespace std;

int main() {
	int rank, size, ierr;
	MPI_Comm comm;

	comm = MPI_COMM_WORLD;

	MPI_Init(NULL, NULL);
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	float T, dt, dx, mm, t, x, y;
	int m, n, m1, n1, n2, n3, n0, pdd, gg, pop, sepp, startt;
	int i, j, k, g,count, fi = 0;
	double t1, t2, t3, t4, t0,qqq;
	double pi = 3.14159265358979323846;
	double aveg[102];
	//MPI_Barrier(MPI_COMM_WORLD);
	t0 = MPI_Wtime();

	for (i = 0; i < 102; i++) {
		aveg[i] = 0.0;
	}

	pdd = 1;
	m = 2304;
	n = 5;
	mm = 2304;
	pop = 2304 % size;
	dx = 2.0 / 2305.0;  // x,y \in (-1,1)
	dt = 0.2 / 2306;
	n0 = 0;//t=0;
	mm = (1.0 / 3.0) / (0.2 / 2306);
	n1 = ceil(mm);   //t=1/3
	mm = (2.0 / 3.0) / (0.2 / 2306);
	n2 = ceil(mm); //t=2/3
	mm = (3.0 / 3.0) / (0.2 / 2306);
	n3 = ceil(mm); //t=3/3


	// how many lines each core needs to caclualte
	int sep;
	sep = (2304) / size;
	float dtx = 0.1 * 0.1;
	if ((2304) % size != 0) { sep++; }


	//vector<vector<float> >    inital2(sep + 2, vector<float>(2306)), inital(sep+2, vector<float>(2306)),calcu(sep+2, vector<float>(2306));
	double inital2[sep + 2][2306], inital[sep + 2][2306], calcu[sep + 2][2306];
	double head[2307], tail[2307];

	m = 2305;
	//sepp is how many lines needs to be calculated
	sepp = sep;
	pop = 2304 % size;
	//initialize for each core, last core is still fine, since just initial
	if (rank + 1 > pop) { sepp = sep - 1; }

	sepp = 23;
	sep = 23;

	//startt is the 0line of the core;
	if (rank + 1 > pop) { startt = rank * sep - (rank)+pop; }
	if (rank + 1 <= pop) { startt = rank * sep; }


	if (rank < size) {


		for (i = 0; i <= sepp + 1; i++) {
			for (j = 1; j <= m - 1; j++) {
				x = -1.0 + dx * j;

				y = -1.0 + dx * (i + startt);
				inital[i][j] = exp(-40 * ((x - 0.4) * (x - 0.4) + y * y));
				//cout << "sep " << inital[i][j] << endl;
			}
			inital[i][0] = 0;
			inital[i][m] = 0;

		}
	}



	if (rank == 0) {
		for (j = 0; j <= m; j++) {
			inital[0][j] = 0.0;
		}
	}

	if (rank == size - 1) {
		for (j = 0; j <= m; j++) {
			inital[sepp][j] = 0.0;
		}
	}
	for (i = 0; i <= sepp + 1; i++) {
		inital[i][0] = 0.0;
		inital[i][m] = 0.0;
	}

	for (i = 0; i <= sepp + 1; i++) {
		for (j = 0; j <= m; j++) {
			inital2[i][j] = inital[i][j];
		}
	}
	// Relation between Actual i >> cacu f(i):  f(i)=rank*sep+i  



	//MPI_Send(&(c1[0]), 7, MPI_INT, 0, 1, MPI_COMM_WORLD);
	//MPI_Recv(&(c2), 7, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  

	if (rank == 0) {
		cout << "i am n1,n2,n3 " << n1 << "  " << n2 << "  " << n3 << "  " << endl;
		cout << "sep " << sep << endl;
		cout << "dtx    dx" << dtx << "  " << dx << endl;
	}



	for (k = 2; k <= n3 + 1; k++) { // n3 steps on time variable
		pdd = 1;
		if ((rank == 0) && (k % 100 == 0)) {
			cout << "k  " << k << endl;
			//cout << "size  " << inital.size() << endl;
		}
		//initalize inital 000
		if (rank == 0) {
			for (j = 0; j <= m; j++) {
				inital[0][j] = 0.0;
			}
		}

		for (i = 0; i <= sepp + 1; i++) {
			inital[i][0] = 0.0;
			inital[i][m] = 0.0;
		}



		//calculation:  A.swap(B);
		for (i = 1; i <= sepp; i++) {
			for (j = 1; j <= m - 1; j++) {
				calcu[i][j] = 2.0 * inital[i][j] - inital2[i][j] + dtx * (inital[i + 1][j] + inital[i - 1][j] + inital[i][j + 1] + inital[i][j - 1] - 4.0 * inital[i][j]);
				
			}
		}
		for (g = 1; g <= 107; g++) {
			for (i = 1; i <= sepp; i++) {
				for (j = 1; j <= m - 1; j++) {
					calcu[i][j] = 2.0 * inital[i][j] - inital2[i][j] + dtx * (inital[i + 1][j] + inital[i - 1][j] + inital[i][j + 1] + inital[i][j - 1] - 4.0 * inital[i][j]);
				}
			}

		}
		


		if (k == n3 + 1) { continue; }

		//make head tail:

		for (i = 1; i < m; i++) {
			head[i] = calcu[1][i];
			tail[i] = calcu[sepp][i];
		}

		//communication:
		//head to his prev
		tail[0] = 0;
		tail[m] = 0;
		head[0] = 0;
		head[m] = 0;


		//renew inital inital2
		tail[0] = 0;
		tail[m] = 0;
		head[0] = 0;
		head[m] = 0;


		for (i = 0; i <= m; i++) {
			calcu[0][i] = head[i];
			calcu[sepp + 1][i] = tail[i];
		}

		if (rank == size - 1) {
			for (i = 0; i <= m; i++) {
				calcu[0][i] = head[i];
				calcu[sepp + 1][i] = 0;
			}
		}

		if (rank == 0) {
			for (i = 0; i <= m; i++) {
				calcu[0][i] = 0;
				calcu[sepp + 1][i] = tail[i];
			}
		}



		//initalize inital 000
		if (rank == 0) {
			for (j = 0; j <= m; j++) {
				inital[0][j] = 0.0;
			}
		}
		if (rank == size - 1) {
			for (j = 0; j <= m; j++) {
				inital[sepp][j] = 0.0;
			}
		}
		for (i = 0; i <= sepp + 1; i++) {
			inital[i][0] = 0.0;
			inital[i][m] = 0.0;
		}


		for(g = 1; g <= 107; g++) {
			for (i = 0; i <= sepp + 1; i++) {
				for (j = 0; j <= m; j++) {
					inital2[i][j] = inital2[i][j]+0.0000000000001;
				}
			}

			for (i = 0; i <= sepp + 1; i++) {
				for (j = 0; j <= m; j++) {
					inital[i][j] = inital[i][j]+0.000000000000001;
				}
			}
		
		}


		
		//initalize inital 000
		if (rank == 0) {
			for (j = 0; j <= m; j++) {
				inital[0][j] = 0.0;
			}
		}
		if (rank == size - 1) {
			for (j = 0; j <= m; j++) {
				inital[sepp][j] = 0.0;
			}
		}
		for (i = 0; i <= sepp + 1; i++) {
			inital[i][0] = 0.0;
			inital[i][m] = 0.0;
		}


	}

	//MPI_Barrier(MPI_COMM_WORLD);
	t4 = MPI_Wtime();

	if (rank == 0) {
		cout << " t0  " << t0 << endl;
		cout << " t1  " << t1 << endl;
		cout << " t2  " << t2 << endl;
		cout << " t3  " << t3 << endl;
		cout << " t4  " << t4 << endl;

	}

	MPI_Finalize();



}





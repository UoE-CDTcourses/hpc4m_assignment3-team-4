#include <iostream>
#include <mpi.h>
#include <math.h>  
#include <vector>
using namespace std;

// We use 10x10 cores, each cores 230x230;

int main(){
  int rank, size, ierr;
  MPI_Comm comm;

  comm  = MPI_COMM_WORLD;
  
  MPI_Init(NULL,NULL);
  MPI_Comm_rank(comm, &rank);            
  MPI_Comm_size(comm, &size);
  
  float T, dt, dx,mm,t,x,y;
  int m, n,m1,n1,n2,n3,n0,pdd,gg,pop,sepp,xstart,ystart;
  int i, j,k, count,fi = 0;
  double t1, t2, t3,t4,t0;
  double pi = 3.14159265358979323846;
  double aveg[105][105];
  //MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  

  if (rank <= 99) {

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
	  int sepx, sepy;
	  sepx = 230;
	  sepy = 230;
	  if ((rank + 1) % 10 == 0) {
		  sepx = 234;
	  }
	  if ((rank) >= 90) {
		  sepy = 234;
	  }

	  float dtx = 0.1 * 0.1;

	  //vector<vector<float> >    inital2(sep + 2, vector<float>(2306)), inital(sep+2, vector<float>(2306)),calcu(sep+2, vector<float>(2306));
	  double inital2[236][236], inital[236][236], calcu[236][236];



	  double xl[236], xr[236], yu[236], yd[236];

	  m = 2305;


	  //n1=50; n2=70; n3=100;

	  //sepp is how many lines needs to be calculated
	  //initialize for each core, last core is still fine, since just initial

	  //xstart  ystart;
	  xstart = (rank % 10) * 230;
	  ystart = (rank / 10) * 230;



	  if (rank < size) {


		  for (i = 0; i <= sepx + 1; i++) {
			  for (j = 0; j <= sepy; j++) {
				  x = -1.0 + dx * (i + xstart);

				  y = -1.0 + dx * (j + ystart);
				  inital[i][j] = exp(-40 * ((x - 0.4) * (x - 0.4) + y * y));
				  //cout << "sep " << inital[i][j] << endl;
			  }

		  }
	  }


	  //initalize inital 000
	  if (rank <= 9) {
		  for (i = 0; i <= sepx + 1; i++) {
			  inital[i][0] = 0.0;
		  }
	  }
	  if (rank >= 90) {
		  for (i = 0; i <= sepx + 1; i++) {
			  inital[i][sepy + 1] = 0.0;
		  }
	  }
	  if (rank % 10 == 0) {
		  for (j = 0; j <= sepy + 1; j++) {
			  inital[0][j] = 0.0;
		  }
	  }
	  if (rank % 10 == 9) {
		  for (j = 0; j <= sepy + 1; j++) {
			  inital[sepx + 1][j] = 0.0;
		  }
	  }





	  for (i = 0; i <= sepx + 1; i++) {
		  for (j = 0; j <= sepy + 1; j++) {
			  inital2[i][j] = inital[i][j];
		  }
	  }
	  // Relation between Actual i >> cacu f(i):  f(i)=rank*sep+i  



	  //MPI_Send(&(c1[0]), 7, MPI_INT, 0, 1, MPI_COMM_WORLD);
	  //MPI_Recv(&(c2), 7, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  

	  if (rank == 0) {
		  cout << "i am n1,n2,n3 " << n1 << "  " << n2 << "  " << n3 << "  " << endl;
	  }


	  for (k = 2; k <= n3 + 1; k++) { // n3 steps on time variable
		  pdd = 1;
		  if ((rank == 0) && (k % 100 == 0)) {
			  cout << "k  " << k << endl;
			  //cout << "size  " << inital.size() << endl;
		  }
		  //initalize inital 000
		  if (rank <= 9) {
			  for (i = 0; i <= sepx + 1; i++) {
				  inital[i][0] = 0.0;
			  }
		  }
		  if (rank >= 90) {
			  for (i = 0; i <= sepx + 1; i++) {
				  inital[i][sepy + 1] = 0.0;
			  }
		  }
		  if (rank % 10 == 0) {
			  for (j = 0; j <= sepy + 1; j++) {
				  inital[0][j] = 0.0;
			  }
		  }
		  if (rank % 10 == 9) {
			  for (j = 0; j <= sepy + 1; j++) {
				  inital[sepx + 1][j] = 0.0;
			  }
		  }



		  //calculation:  A.swap(B);
		  for (i = 1; i <= sepx; i++) {
			  for (j = 1; j <= sepy; j++) {
				  calcu[i][j] = 2.0 * inital[i][j] - inital2[i][j] + dtx * (inital[i + 1][j] + inital[i - 1][j] + inital[i][j + 1] + inital[i][j - 1] - 4.0 * inital[i][j]);
			  }
		  }


		  if (k == n3 + 1) { continue; }

		  //make xl xr yu yd:

		  for (i = 1; i < sepx; i++) {
			  yd[i] = calcu[i][1];
			  yu[i] = calcu[i][sepy];
		  }
		  for (j = 1; j < sepy; j++) {
			  xl[j] = calcu[1][j];
			  xr[j] = calcu[sepx][j];
		  }



		  //communication:
		  //head to his prev
		  xl[0] = 0;
		  xl[235] = 0;
		  xr[0] = 0;
		  xr[235] = 0;
		  yu[0] = 0;
		  yu[235] = 0;
		  yd[0] = 0;
		  yd[235] = 0;

		  //send yd
		  if ((rank >= 10)) {
			  MPI_Send(&(yd[0]), 236, MPI_DOUBLE, (rank - 10), k, MPI_COMM_WORLD);
		  }

		  //send yu
		  if ((rank < 90)) {
			  MPI_Send(&(yu[0]), 236, MPI_DOUBLE, (rank + 10), k, MPI_COMM_WORLD);
		  }

		  //send xl
		  if ((rank % 10 > 0)) {
			  MPI_Send(&(xl[0]), 236, MPI_DOUBLE, (rank - 1), k, MPI_COMM_WORLD);
		  }
		  //send xr
		  if ((rank % 10 < 9)) {
			  MPI_Send(&(xr[0]), 236, MPI_DOUBLE, (rank + 1), k, MPI_COMM_WORLD);
		  }

		  //MPI_Recv(&pdd, 1, MPI_INT, rank - 1, -1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		  //receive his head
		  //receive yd
		  if ((rank >= 10)) {
			  MPI_Recv(&(yd[0]), 236, MPI_DOUBLE, (rank - 10), k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		  }

		  //rec yu
		  if ((rank < 90)) {
			  MPI_Recv(&(yu[0]), 236, MPI_DOUBLE, (rank + 10), k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		  }

		  //rec xl
		  if ((rank % 10 > 0)) {
			  MPI_Recv(&(xl[0]), 236, MPI_DOUBLE, (rank - 1), k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		  }
		  //rec xr
		  if ((rank % 10 < 9)) {
			  MPI_Recv(&(xr[0]), 236, MPI_DOUBLE, (rank + 1), k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		  }




		  //renew inital inital2
		  xl[0] = 0;
		  xl[235] = 0;
		  xr[0] = 0;
		  xr[235] = 0;
		  yu[0] = 0;
		  yu[235] = 0;
		  yd[0] = 0;
		  yd[235] = 0;

		  for (i = 1; i < sepx; i++) {
			  calcu[i][0] = yd[i];
			  calcu[i][sepy + 1] = yu[i];
		  }
		  for (j = 1; j < sepy; j++) {
			  calcu[0][j] = xl[j];
			  calcu[sepx + 1][j] = xr[j];
		  }






		  //initalize inital 000
		  if (rank <= 9) {
			  for (i = 0; i <= sepx + 1; i++) {
				  inital[i][0] = 0.0;
			  }
		  }
		  if (rank >= 90) {
			  for (i = 0; i <= sepx + 1; i++) {
				  inital[i][sepy + 1] = 0.0;
			  }
		  }
		  if (rank % 10 == 0) {
			  for (j = 0; j <= sepy + 1; j++) {
				  inital[0][j] = 0.0;
			  }
		  }
		  if (rank % 10 == 9) {
			  for (j = 0; j <= sepy + 1; j++) {
				  inital[sepx + 1][j] = 0.0;
			  }
		  }

		  //renew inital inital2
		  for (i = 0; i <= sepx + 1; i++) {
			  for (j = 0; j <= sepy + 1; j++) {
				  inital2[i][j] = inital[i][j];
			  }
		  }
		  for (i = 0; i <= sepx + 1; i++) {
			  for (j = 0; j <= sepy + 1; j++) {
				  inital[i][j] = calcu[i][j];
			  }
		  }


		  //initalize inital 000
		  if (rank <= 9) {
			  for (i = 0; i <= sepx + 1; i++) {
				  inital[i][0] = 0.0;
			  }
		  }
		  if (rank >= 90) {
			  for (i = 0; i <= sepx + 1; i++) {
				  inital[i][sepy + 1] = 0.0;
			  }
		  }
		  if (rank % 10 == 0) {
			  for (j = 0; j <= sepy + 1; j++) {
				  inital[0][j] = 0.0;
			  }
		  }
		  if (rank % 10 == 9) {
			  for (j = 0; j <= sepy + 1; j++) {
				  inital[sepx + 1][j] = 0.0;
			  }
		  }



		  //Output
		  if ((k == n1) || (k == n2) || (k == n3)) {
			  // Conclude results. Send to core 0, mark as true row value

			  if (k == n1) {
				  //MPI_Barrier(MPI_COMM_WORLD);
				  t1 = MPI_Wtime();
			  }
			  if (k == n2) {
				  //MPI_Barrier(MPI_COMM_WORLD);
				  t2 = MPI_Wtime();
			  }

			  for (i = 0; i <= 100; i++) {
				  //cout << "real i   k  " << i << "  " << k << endl;
				  for (j = 0; j <= 100; j++) {
					  aveg[i][j] = aveg[i][j] = 0.00000001;

				  }
			  }
			  //sep=22, sepp=21 xor 22
			  cout << "rank k  " << rank << "  " << k << endl;
			  if (rank == 0) {
				  for (i = 1; i <= 229; i++) {
					  //cout << "real i   k  " << i << "  " << k << endl;
					  for (j = 1; j <= 229; j++) {
						  aveg[i / 23][j / 23] += calcu[j][i];
						  //cout << calcu[i][j] << " ";
					  }
				  }

				  for (i = 0; i <= 99; i++) {
					  //cout << "real i   k  " << i << "  " << k << endl;
					  for (j = 0; j <= 99; j++) {
						  aveg[i][j] = aveg[i][j] / 100.0;
						  cout << calcu[i][j] << " ";
					  }
				  }
				  cout << endl;
				  pdd = 0;
				  MPI_Send(&(pdd), 1, MPI_INT, 1, -1, MPI_COMM_WORLD);
				  pdd = 1;
			  }
			  if (rank > 0) {
				  MPI_Recv(&pdd, 1, MPI_INT, rank - 1, -1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				  if (pdd == 0) {
					  for (i = 1; i <= 229; i++) {
						  //cout << "real i   k  " << i << "  " << k << endl;
						  for (j = 1; j <= 229; j++) {
							  aveg[i / 23][j / 23] += calcu[j][i];
							  //cout << calcu[i][j] << " ";
						  }
					  }

					  for (i = 0; i <= 99; i++) {
						  //cout << "real i   k  " << i << "  " << k << endl;
						  for (j = 0; j <= 99; j++) {
							  aveg[i][j] = aveg[i][j] / 100.0;
							  cout << calcu[i][j] << " ";
						  }
					  }
					  cout << endl;
					  pdd = 0;
					  if (rank <= 98) {
						  MPI_Send(&(pdd), 1, MPI_INT, rank + 1, -1, MPI_COMM_WORLD);
					  }

					  pdd = 1;

				  }
			  }


		  }

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

 



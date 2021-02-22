#include <iostream>
#include <mpi.h>
#include <math.h>  
#include <vector>
using namespace std;

int main(){
  int rank, size, ierr;
  MPI_Comm comm;

  comm  = MPI_COMM_WORLD;
  
  MPI_Init(NULL,NULL);
  MPI_Comm_rank(comm, &rank);            
  MPI_Comm_size(comm, &size);
  
  float T, dt, dx,mm,t,x,y;
  int m, n,m1,n1,n2,n3,n0,pdd,gg,pop,sepp,startt;
  int i, j,k, count,fi = 0;
  double t1, t2, t3,t4,t0;
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
  dt = 0.2/2306;
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
  float dtx = 0.1*0.1;
  if ((2304 ) % size != 0) { sep++; }


  //vector<vector<float> >    inital2(sep + 2, vector<float>(2306)), inital(sep+2, vector<float>(2306)),calcu(sep+2, vector<float>(2306));
  double inital2[2306][sep + 2], inital[2306][sep + 2], calcu[2306][sep + 2];



  double head[2307], tail[2307];

  m = 2305;


  //n1=50; n2=70; n3=100;

  //sepp is how many lines needs to be calculated
  sepp = sep;
  pop = 2304 % size;
  //initialize for each core, last core is still fine, since just initial
  if (rank + 1 > pop) {   sepp = sep - 1;  }

  //startt is the 0line of the core;
  if (rank + 1 > pop) { startt=rank*sep-(rank)+pop; }
  if (rank + 1 <= pop) { startt=rank*sep ; }


  if (rank < size ) {


	  for (i = 0; i <= sepp + 1; i++) {
		  for (j = 1; j <= m-1; j++) {
			  x = -1.0+dx * j;

			  y = -1.0+dx * (i + startt);
			  inital[j][i] = exp(-40 * ((x - 0.4) * (x - 0.4) + y * y));
			  //cout << "sep " << inital[i][j] << endl;
		  }
		  inital[0][i] = 0;
		  inital[m][i] = 0;

	  }
  }

 
  //initalize inital 000
  if (rank == 0) {
	  for (j = 0; j <= m; j++) {
		  inital[j][0] = 0.0;
	  }
  }

  if (rank == size - 1) {
	  for (j = 0; j <= m; j++) {
		  inital[j][sepp] = 0.0;
	  }
  }
  for (i = 0; i <= sepp + 1; i++) {
	  inital[0][i] = 0.0;
	  inital[m][i] = 0.0;
  }

  for (i = 0; i <= sepp + 1; i++) {
	  for (j = 0; j <= m ; j++) {
		  inital2[j][i] = inital[j][i];
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
	  if ((rank == 0)&&(k%100==0)) {
		  cout << "k  " << k << endl;
		  //cout << "size  " << inital.size() << endl;
	  }
	  //initalize inital 000
	  if (rank == 0) {
		  for (j = 0; j <= m; j++) {
			  inital[j][0] = 0.0;
		  }
	  }

	  if (rank == size - 1) {
		  for (j = 0; j <= m; j++) {
			  inital[j][sepp] = 0.0;
		  }
	  }
	  for (i = 0; i <= sepp + 1; i++) {
		  inital[0][i] = 0.0;
		  inital[m][i] = 0.0;
	  }



	  //calculation:  A.swap(B);
	  for (i = 1; i <= sepp; i++) {
		  for (j = 1; j <= m - 1; j++) {
			  calcu[j][i] = 2.0 * inital[j][i] - inital2[j][i] + dtx * (inital[j + 1][i] + inital[j - 1][i] + inital[j][i + 1] + inital[j][i - 1] - 4.0 * inital[j][i]);
		  }
	  }


	  if (k == n3+1) { continue; }

	  //make head tail:

	  for (i = 1; i < m; i++) {
		  head[i] = calcu[i][1];
		  tail[i] = calcu[i][sepp];
	  }

	  //communication:
	  //head to his prev
	  tail[0] = 0;
	  tail[m] = 0;
	  head[0] = 0;
	  head[m] = 0;
	  if (rank > 0) {
		  MPI_Send(&(head[0]), m +1, MPI_DOUBLE, (rank - 1), k, MPI_COMM_WORLD);
	  }

	  //tail to his back
	  if (rank < size - 1) {
		  MPI_Send(&(tail[0]), m + 1, MPI_DOUBLE, (rank + 1), k, MPI_COMM_WORLD);
	  }

	  //receive his head
	  if (rank > 0) {
		  MPI_Recv(&(head[0]), m + 1, MPI_DOUBLE, rank - 1, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  }

	  if (rank < size - 1) {
		  MPI_Recv(&(tail[0]), m + 1, MPI_DOUBLE, rank + 1, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  }

	  //renew inital inital2
	  tail[0] = 0;
	  tail[m] = 0;
	  head[0] = 0;
	  head[m] = 0;


	  for (i = 0; i <= m; i++) {
		  calcu[i][0] = head[i];
		  calcu[i][sepp+1] = tail[i];
	  }

	  if (rank == size - 1) {
		  for (i = 0; i <= m; i++) {
			  calcu[i][0] = head[i];
			  calcu[i][sepp+1] = 0;
		  }
	  }

	  if (rank == 0) {
		  for (i = 0; i <= m; i++) {
			  calcu[i][0] =0;
			  calcu[i][sepp+1] = tail[i];
		  }
	  }




	  //initalize inital 000
	  if (rank == 0) {
		  for (j = 0; j <= m; j++) {
			  inital[j][0] = 0.0;
		  }
	  }

	  if (rank == size - 1) {
		  for (j = 0; j <= m; j++) {
			  inital[j][sepp] = 0.0;
		  }
	  }
	  for (i = 0; i <= sepp + 1; i++) {
		  inital[0][i] = 0.0;
		  inital[m][i] = 0.0;
	  }




	  for (i = 0; i <= sepp + 1; i++) {
		  for (j = 0; j <= m; j++) {
			  inital2[j][i] = inital[j][i];
		  }
	  }

	  for (i = 0; i <= sepp + 1; i++) {
		  for (j = 0; j <= m; j++) {
			  inital[j][i] = calcu[j][i];
		  }
	  }
	  //initalize inital 000
	  if (rank == 0) {
		  for (j = 0; j <= m; j++) {
			  inital[j][0] = 0.0;
		  }
	  }

	  if (rank == size - 1) {
		  for (j = 0; j <= m; j++) {
			  inital[j][sepp] = 0.0;
		  }
	  }
	  for (i = 0; i <= sepp + 1; i++) {
		  inital[0][i] = 0.0;
		  inital[m][i] = 0.0;
	  }




	  if ((k == n1) || (k == n2)  || (k == n3)) { // Choose which to take
		  // Conclude results. Send to core 0, mark as true row value

		  if (k == n1) {
			  //MPI_Barrier(MPI_COMM_WORLD);
		      t1 = MPI_Wtime();
		  }
		  if (k == n2) {
			  //MPI_Barrier(MPI_COMM_WORLD);
			  t2 = MPI_Wtime();
		  }

		  for (i = 0; i < 102; i++) {
			  aveg[i] = 0.0000000001;
		  }
		  //sep=22, sepp=21 xor 22
		  if (rank == 0) {
			  for (i = 1; i <= sepp; i++) {
				  //cout << "real i   k  " << i << "  " << k << endl;
				  for (j = 0; j <= m; j++) {
					  gg = j / 23;
					  aveg[gg] += calcu[j][i];
					  //cout << calcu[i][j] << " ";
				  }
			  }
			  cout << "rank k  " << rank<< "  " <<k<<  endl;
			  for (j = 0; j <= 100; j++) {
				  aveg[j] = aveg[j] / (23.0 * sepp);
				  cout << aveg[j] << " ";
				  //cout << calcu[10][j] << " ";
			  }
			  cout << endl;
			  pdd = 0;
			  MPI_Send(&(pdd), 1, MPI_INT, 1, -1, MPI_COMM_WORLD);
			  pdd = 1;
		  }

		  if (rank > 0 && rank < size - 1) {
			  MPI_Recv(&pdd, 1, MPI_INT, rank - 1, -1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			  //cout << "xxxxx   k  " << rank << "  " << k << endl;
			  if (pdd == 0) {
				  for (i = 1; i <= sepp; i++) {
					  //cout << "real i   k  " << i+rank*sep << "  " << k << endl;
					  for (j = 0; j <= m; j++) {
						  gg = j / 23;
						  aveg[gg] += calcu[j][i];
						  //cout << calcu[i][j] << " ";
					  }
				  }
				  cout << "rank k  " << rank << "  " << k << endl;
				  for (j = 0; j <= 100; j++) {
					  aveg[j] = aveg[j] / (23.0 * sepp);
					  cout << aveg[j] << " ";
					  //cout << calcu[10][j] << " ";
				  }
				  cout << endl;




				  pdd = 0;
				  MPI_Send(&(pdd), 1, MPI_INT, rank+1, -1, MPI_COMM_WORLD);
				  pdd = 1;
			  }
		  }
		  if (rank == size - 1) {
			  MPI_Recv(&pdd, 1, MPI_INT, rank - 1, -1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			  cout << "xxxxx   k  " << rank << "  " << k << endl;
			  if (pdd == 0) {
				  for (i = 1; i <= sepp; i++) {
					  //cout << "real i   k  " << i + rank * sep << "  " << k << endl;
					  for (j = 0; j <= m; j++) {
						  gg = j / 23;
						  aveg[gg] += calcu[j][i];
						  //cout << calcu[i][j] << " ";
					  }
				  }
				  for (j = 0; j <= 100; j++) {
					  aveg[j] = aveg[j] / (23.0 * sepp);
					  cout << aveg[j] << " ";
					  //cout << calcu[10][j] << " ";
				  }
				  cout << endl;
				  pdd = 1;

			  }
		  }



		  int coree;

			  /*
			  for (i = sep + 1; i <= m - 1; i++) {
				  coree = i / sep;
				  //MPI_Recv(&(head[1]), m - 1, MPI_FLOAT, coree, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				  head[0] = 0;
				  head[m] = 0;

				  for (j = 0; j <= m; j++) {
					  //cout << head[j] << " ";
				  }
				  cout << "xxxxx  " << i << "  " << endl;
				  cout << " " << endl;


			  }*/

		
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

 



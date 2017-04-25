#include <math.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
    
//#include "timer.h"

#define N 5
#define P_SIZE ((int) (pow(12, N+1)-1)/11)

/* Barriers */
long double E_1_l = 2.951646;
long double E_4_l = -3590.207247+3592.956166;
long double E_1_4 = -3589.95316882+3592.97102056;
long double E_3_4 = -3589.37845071+3592.97101757;
long double E_c =  2.781326; //2.492666; //3.067482; 

/* Boltzmann constant */
long double kb = 8.6173324e-5;

/* Vibrational frequency and temperature*/
long double v = 1.6e13;
long double T = 1300;
/* T = 1800; */

/* Nearest neighbour distance */
long double nnd =  3.08472680894400e-1;

/* Silicon atom positions in unit cell */
long double spos_Si[6][3] = {
                              {0.33333333,  0.66666667,  0.93750000-1},
                              {0.00000, 0.00000, 0.1875},
                              {0.66666667,  0.33333333,  0.43750000},
                              {0.000000,  0.000000,  0.68750000},
                              {0.33333333,  0.66666667,  0.93750000},
                              {0.00000, 0.00000, 0.1875+1}
                            };

/* lattice scaling parameters and unit cell */
long double a = 1.0104;
long double c = 1.0074;

/* Define two rotation matrix to capture all possible in-plane and out of
   plane transitions */
long double rotm1[3][3] = {
                            {0.5, -0.866025403784439, 0.0},
                            {0.866025403784439, 0.5, 0.0},
                            {0.0, 0.0, 1.0}
                          };


long double rotm2[3][3] = {
                            {-0.5, -0.866025403784439, 0.0},
                            {0.866025403784439, -0.5, 0.0},
                            {0.0, 0.0, 1.0}
                          };

long double cell2[3][3];

long double ipvec[6][3];

long double opvec1[4][3][3];

long double opvec2[4][3][3];

// mat0 is product, mat1, mat2 are things to be multiplied
void mat_mul(long double prod[3][3], long double mat1[3][3], long double mat2[3][3]) {
  int i, j, k;

  // #pragma omp parallel shared(prod, mat1, mat2) private(i, j, k)
  // {
  //   #pragma omp for schedule(static)
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        prod[i][j] = 0.0;
        for (k = 0; k < 3; k++) {
          prod[i][j] += mat1[i][k] * mat2[k][j];
        }
      }
    }
  // }

  return; 
}

// mat0 stores the product, mat2 is the thing to be exponentiated
void mat_pow(long double mat0[3][3], long double mat2[3][3], int power) {
  long double mat1[3][3];
  for (int i = 0; i < 3; i++) {
    memcpy(&mat1[i], &mat2[i], sizeof(mat2[0]));
  }

  // if power == 0, set mat0 to identity matrix
  if (power == 0) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
  	    if (i == j) {
  	      mat0[i][j] = 1.0;
  	    }
    	  else {
    	    mat0[i][j] = 0.0;
    	  }
      }
    }
  }
  
  // if power == 1, set mat0 to mat2
  else if (power == 1) {
    for (int i = 0; i < 3; i++) {
      memcpy(&mat0[i], &mat2[i], sizeof(mat2[0]));
	  }
  }

  for (int p = 1; p < power; p++) {
    // square mat2 and let mat0 equal mat2**2
    mat_mul(mat0, mat1, mat2);
    // let mat1 = mat2**2
    for (int i = 0; i < 3; i++) {
      memcpy(&mat1[i], &mat0[i], sizeof(mat0[0]));
	  }
  }

  return; 
}

void mat_vec_mul(long double prod[3], long double vec[3], long double mat[3][3]) {
  int i, k;
  
  // #pragma omp parallel shared(prod, vec, mat) private (i, k)
  // {
  //   #pragma omp for schedule(static)
    for (i = 0; i < 3; i++) {
      prod[i] = 0.0;
      for (int k = 0; k < 3; k++) {
        prod[i] += vec[k] * mat[k][i];
      }
    }
  // }
  
  return; 
}

/* Run simulation */
// A[i][j] keeps track of allowable transitions between i and j
void BFS(long double A[P_SIZE][5], long double P[P_SIZE][6], int world_rank, int world_size,int m) {
  // omp_set_num_threads(91);
  int rank;
  long int lower, upper, i, j,k,size,width,start,end;
  long double holder[6];
  start = 0;
  end = 1;
   
  for (k = 0;k<m;k++) 
    {
      start += pow(12,k);
      end += pow(12,k+1);
    }
  width = end-start;
  size = world_size < width ? world_size : width;
  k = world_rank > (width-1) ? (long int) world_size: (long int) world_rank;
  for (i = start+(long int) round((k)*(width)/(1.0*size)); i < start+ (long int) round((world_rank+1.0)*width/(1.0*size)); i++) {
    // lower bound of particles we can move to
    lower = ((12*i+1) < P_SIZE) ? (12*i+1) : P_SIZE;
    //  if (lower == P_SIZE) {
    //   break;
    // }
    // upper bound of 12 particles we can move to
    upper = ((12*i+13) <= P_SIZE) ? (12*i+13) : P_SIZE;
    // iterate through the 12 possible lattice movements
    for (j = lower; j < upper; j++) {
      if (P[i][5] == 1.0) {
	// [0] is x position
	P[j][0] = A[j][0]+P[i][0];
	// [5] keeps track of where we visited
	P[j][5] = 1.0;
	//   *ii = *ii+1;
	if ((j % 12) == 0) {
	  P[i][5] += 1.0;
	}
	// [1] is y position
	P[j][1] = A[j][1]+P[i][1];
	// [2] is z position
	P[j][2] = A[j][2]+P[i][2];
	// [3] is time
	P[j][3] = A[j][3]+P[i][3];
	// [4] is probability
	P[j][4] = A[j][4]*P[i][4];
      }
    }
  }
  for (i = start; i < start+width; i++) {
    // lower bound of particles we can move to
    lower = ((12*i+1) < P_SIZE) ? (12*i+1) : P_SIZE;
    //  if (lower == P_SIZE) {
    //   break;
    // }
    // upper bound of 12 particles we can move to
    upper = ((12*i+13) <= P_SIZE) ? (12*i+13) : P_SIZE;
    rank = (int) floor((i-start)*1.0*size/(1.0*width));
    {
      // iterate through the 12 possible lattice movements
      for (j = lower; j < upper; j++) {
	holder[0] = P[j][0];
	holder[1] = P[j][1];
	holder[2] = P[j][2];
	holder[3] = P[j][3];
	holder[4] = P[j][4];
	holder[5] = P[j][5];
	MPI_Bcast(holder,6,MPI_LONG_DOUBLE,rank,MPI_COMM_WORLD);
	P[j][0] = holder[0];
	P[j][1] = holder[1];
	P[j][2] = holder[2];
	P[j][3] = holder[3];
	P[j][4] = holder[4];
	P[j][5] = holder[5];
      }
    }
  }
}
 
int main(int argc, char** argv) {
  
  MPI_Init(NULL, NULL);
  time_t start = time(NULL);
  time_t end = time(NULL);
  long double r1 = v*exp(-E_1_l/kb/T);
  long double r2 = v*exp(-E_4_l/kb/T);
  long double r3 = v*exp(-E_1_4/kb/T);
  long double r4 = v*exp(-E_3_4/kb/T);
  long double r5 = v*exp(-E_c/kb/T);
  long double Q1, Q2;
  long int i, j, k,kk,l,m, width,size,endd,startt;

  long double (*P)[6];
  long double (*P2)[5];
  long double (*A)[5];
  long double rotmf[3][3];
  long double mf[3][3];
  long double vec1[3];
  long double vec2[3];
  long int n_index;
  long double sum, gsum;
  //StartTimer();
  // transition rates
  // starting from cubic site
  Q1 = 6.*r1+3.*r3+3.*r4+3.*r5;
  // starting from hexagonal site
  Q2 = 6.*r2+3.*r3+3.*r4+3.*r5;
  
  P = malloc(P_SIZE * sizeof(long double[6]));
  P2 = malloc(P_SIZE * sizeof(long double[5]));
  A = malloc(P_SIZE * sizeof(long double[5]));
  int *checker;
  int *swtcher;
  
  checker = malloc(P_SIZE * sizeof(int));
  swtcher = malloc(P_SIZE * sizeof(int));

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  
  
  // can't individually place after initialization
  // lattice vector
  long double cell2_duplicate[3][3] = {
            {nnd, 0*a, 0*a},
            {-nnd/2,nnd/2*sqrt(3), 0*a},
            {0*1.0, 0*1.0, 10.086*c*nnd/3.078*a/2.57218587467527*2.51866888630220}
          };
  // in plane vector
  // out of plane vector 1 is down
  // out of plane vector 2 is up
  // 6 for 60 degrees
  while((end-start)<10)
    {
  for (i = 0; i < 6; i++) {
    for (j = 0; j < 3 ; j++) {
	    vec1[j] = cell2_duplicate[0][j];
	  }
    // doing in plane rotation
    mat_pow(rotmf, rotm1, i);
    // multiply by rot matrix to get rotated lattice vector
    mat_vec_mul(vec2, vec1, rotmf);
    for (j = 0; j < 3; j++) {
      ipvec[i][j] = vec2[j];
    }
  }
  // o.o.p transitions
  for (k = 0; k < 4; k++) {
    for (i = 0; i < 3; i++) {

      // remove duplicated vec1
      for (j = 0; j < 3; j++) {
        vec1[j] = spos_Si[k][j] - spos_Si[k+1][j];
      }
  	  mat_pow(rotmf, rotm2, i);
  	  mat_mul(mf, cell2_duplicate, rotmf);
  	  mat_vec_mul(vec2, vec1, mf);
	    for (j = 0; j < 3; j++) {
	      opvec1[k][i][j] = vec2[j];
	    }
    }
  }
  for (k = 0; k < 4; k++) {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        vec1[j] = spos_Si[k+2][j] - spos_Si[k+1][j];
      }
  	  mat_pow(rotmf, rotm2, i);
  	  mat_mul(mf, cell2_duplicate, rotmf);
  	  mat_vec_mul(vec2, vec1, mf);
  	  for (j = 0; j < 3; j++) {
        opvec2[k][i][j] = vec2[j];
      }
    }
  }

  for (i = 0; i < 3; i++) {
    memcpy(&cell2[i], &cell2_duplicate[i], sizeof(cell2_duplicate[0]));
  }

  // initializes all elements to 0
  for (i = 0; i < P_SIZE; i++) {
    // prevents duplicating probabilities
    checker[i] = 0;
    }
  swtcher[0] = 0;
  for (m = 1; m < N+1; m++) {
    startt = 0;
    endd = 1;
   
    for (k = 0;k<m;k++) 
      {
	startt += pow(12,k);
	endd += pow(12,k+1);
      }
    width = endd-startt;
    size = world_size < width ? world_size : width;
    kk = world_rank > (width-1) ? (long int) world_size: (long int) world_rank;
    for (i = startt+(long int) round((kk)*(width)/(1.0*size)); i < startt+(long int) round((world_rank+1.0)*width/(1.0*size)); i++) {
      // keeps track of what level you're at
      swtcher[i] = 0;
    }
  }
  
  //omp_set_num_threads(91);
  for (k = 0; k < 12; k++) {
    n_index = k+1;
    if (n_index >= P_SIZE) {
      i = P_SIZE;
      break;
    }
    // in plane
    if (k < 6) {
      for (j = 0; j < 3; j++) {
	A[n_index][j] = ipvec[k][j];
      }
      A[n_index][4] = r1/Q1;
      A[n_index][3] = 1./Q1;
    }
    // below plane
    else if (k < 9) {
      for (j = 0; j < 3; j++) {
	A[n_index][j] = opvec1[swtcher[i]][k-6][j];
      }
      A[n_index][4] = r3/Q1;
      A[n_index][3] = 1./Q1;
      swtcher[n_index] = (swtcher[i]+3)%4;
    }
    // above plane
    else {
      for (j = 0; j < 3; j++) {
	A[n_index][j] = opvec2[swtcher[i]][k-9][j];
      }
      A[n_index][4] = r4/Q1;
      A[n_index][3] = 1./Q1;
      swtcher[n_index] = (swtcher[i]+1)%4;
    }
  }
   
  for (m = 1; m < N; m++) {
    startt = 0;
    endd = 1;
   
    for (k = 0;k<m;k++) 
      {
	startt += pow(12,k);
	endd += pow(12,k+1);
      }
    width = endd-startt;
    size = world_size < width ? world_size : width;
    kk = world_rank > (width-1) ? (long int) world_size: (long int) world_rank;
    for (i = startt+(long int) round((kk)*(width)/(1.0*size)); i < startt+ (long int) round((world_rank+1.0)*width/(1.0*size)); i++) {
      if (swtcher[i] == 0 || swtcher[i] == 2) {
        for (k = 0; k < 12; k++) {
  	      n_index = 12*i+k+1;
  	      if (n_index >= P_SIZE) {
      		  i = P_SIZE;
      		  break;
      		}
          // in plane
  	      if (k < 6) {
      		  for (j = 0; j < 3; j++) {
  		        A[n_index][j] = ipvec[k][j];
  		      }
      		  A[n_index][4] = r1/Q1;
      		  A[n_index][3] = 1./Q1;
  		    }
          // below plane
  	      else if (k < 9) {
      		  for (j = 0; j < 3; j++) {
    		      A[n_index][j] = opvec1[swtcher[i]][k-6][j];
    		    }
      		  A[n_index][4] = r3/Q1;
      		  A[n_index][3] = 1./Q1;
      		  swtcher[n_index] = (swtcher[i]+3)%4;
  		    }
          // above plane
  	      else {
      		  for (j = 0; j < 3; j++) {
  	  	      A[n_index][j] = opvec2[swtcher[i]][k-9][j];
  		      }
      		  A[n_index][4] = r4/Q1;
      		  A[n_index][3] = 1./Q1;
      		  swtcher[n_index] = (swtcher[i]+1)%4;
  	    	}
  	    }
      }
      else {
        for (k = 0; k < 12; k++) {
  	      n_index = 12*i+k+1;
          if (n_index >= P_SIZE) {
      		  i = P_SIZE;
  	        break;
  	      } 
          if (k < 6) {
        	  for (j = 0; j < 3; j++) {
  	          A[n_index][j] = ipvec[k][j];
  	        }
      		  A[n_index][4] = r2/Q2;
      		  A[n_index][3] = 1./Q2;
      		}
          else if (k < 9) {
  	        for (j = 0; j < 3; j++) {
    		      A[n_index][j] = opvec1[swtcher[i]][k-6][j];
    		    }
      		  A[n_index][4] = r4/Q2;
      		  A[n_index][3] = 1./Q2;
      		  swtcher[n_index] = (swtcher[i]+3)%4;
  	      }
          else {
        	  for (j = 0; j < 3; j++) {
    		      A[n_index][j] = opvec2[swtcher[i]][k-9][j];
    		    }
      		  A[n_index][4] = r3/Q2;
      		  A[n_index][3] = 1./Q2;
      		  swtcher[n_index] = (swtcher[i]+1)%4;
  	      }
        }
      }
    }
  }
  P[0][1] = 0.0;
  P[0][2] = 0.0;
  P[0][3] = 0.0;
  P[0][4] = 1.0;
  P[0][5] = 1.0;
  //ii = 0;
  
  if (world_rank == 0){
    printf("%Lf\n", P[0][4]);
  }
  for (i = 0; i < N; i++) {
  BFS(A, P,world_rank,world_size,i);
  }
  
  //printf("%ld\n",ii);
  n_index = N % 2;
  
  for (i = 0; i < P_SIZE; i++) {
    P2[i][0] = P[i][0];
    P2[i][1] = P[i][1];
    P2[i][2] = P[i][2];
    P2[i][3] = P[i][3];
    P2[i][4] = 0.0;
    sum = 0.0;

    //omp_set_num_threads(91);
    width = P_SIZE - i;
    size = world_size < width ? world_size : width;
    kk = world_rank > (width-1) ? (long int) world_size: (long int) world_rank;
    for (j = i+(long int) round((kk)*(width)/(1.0*size)); j < i+ (long int) round((world_rank+1.0)*width/(1.0*size)); j++) {
      if (((checker[j] == 0 && fabsl(P[j][0] - P[i][0])<0.0000001) &&
	   (fabsl(P[j][1] - P[i][1])<0.0000001 && fabsl(P[j][2] - P[i][2])<0.0000001)) &&
          (fabsl(P[j][3] - P[i][3])<0.0000001)) {
  	      sum += P[j][4];
  	      checker[j] = 1;
  	  }
    }
    MPI_Reduce(&sum,&gsum,1,MPI_LONG_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); 
    P2[i][4] = gsum;
  }
  
  //if (P2[i][4] > 0.0) {
  //	  printf("%Lf %Lf %Lf %Lf %.11Lf\n", P2[i][0], P2[i][1], P2[i][2], P2[i][3], P2[i][4]);
  //	}

  end = time(NULL);
  l++;
    }
  if (world_rank == 0) {
    printf("\n%.11Lf %ld\n", (long double)((end-start)/(1.0*l)*1000000.),l);
  }
  free(P);
  free(P2);
  free(A);
  free(checker);
  free(swtcher);
  // Finalize the MPI environment.
  MPI_Finalize();

  return 0;
}

#include <math.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N 1
#define P_SIZE ((int) (pow(12, N+1)-1)/11)

/* Barriers */
long double E_1_l = 2.951646;
long double E_4_l = -3590.207247+3592.956166;
long double E_1_4 = -3589.95316882+3592.97102056;
long double E_3_4 = -3589.37845071+3592.97101757;
long double E_c = 3.067482; 
long double r1,r2,r3,r4,r5;
long double Q1, Q2;


/* Boltzmann constant */
long double kb = 8.6173324e-5;

/* Vibrational frequency and temperature*/
long double v = 1.6e13;
long double T = 1300;
/* T = 1800; */

/* transition rates */
long double r1, r2, r3, r4, r5;
long double Q1, Q2;

/* Nearest neighbour distance */
long double nnd = 3.0847268089e-1;

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

  prod[0][0] = mat1[0][0]*mat2[0][0] + mat1[0][1]*mat2[1][0] + mat1[0][2]*mat2[2][0];
  prod[0][1] = mat1[0][0]*mat2[0][1] + mat1[0][1]*mat2[1][1] + mat1[0][2]*mat2[2][1];
  prod[0][2] = mat1[0][0]*mat2[0][2] + mat1[0][1]*mat2[1][2] + mat1[0][2]*mat2[2][2];
  prod[1][0] = mat1[1][0]*mat2[0][0] + mat1[1][1]*mat2[1][0] + mat1[1][2]*mat2[2][0];
  prod[1][1] = mat1[1][0]*mat2[0][1] + mat1[1][1]*mat2[1][1] + mat1[1][2]*mat2[2][1];
  prod[1][2] = mat1[1][0]*mat2[0][2] + mat1[1][1]*mat2[1][2] + mat1[1][2]*mat2[2][2];
  prod[2][0] = mat1[2][0]*mat2[0][0] + mat1[2][1]*mat2[1][0] + mat1[2][2]*mat2[2][0];
  prod[2][1] = mat1[2][0]*mat2[0][1] + mat1[2][1]*mat2[1][1] + mat1[2][2]*mat2[2][1];
  prod[2][2] = mat1[2][0]*mat2[0][2] + mat1[2][1]*mat2[1][2] + mat1[2][2]*mat2[2][2];

  return; 
}

void mat_vec_mul(long double prod[3], long double vec[3], long double mat[3][3]) { 
 
  prod[0] = mat[0][0]*vec[0] + mat[0][1]*vec[1] + mat[0][2]*vec[2];
  prod[1] = mat[1][0]*vec[0] + mat[1][1]*vec[1] + mat[1][2]*vec[2];
  prod[2] = mat[2][0]*vec[0] + mat[2][1]*vec[1] + mat[2][2]*vec[2];

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

/* Run simulation */
void BFS(long double P[P_SIZE][5], long int index, int swtch, long int *index2) {

  *index2 = *index2+1;
  if (12*index+1 >= P_SIZE) {
    return;
  }

  int i, ij;
  // new index
  int n_index;

  if (swtch == 0 || swtch == 2) {
    // #pragma omp parallel shared(vec1, P) private(ij)
    // {
    //   #pragma omp parallel for schedule(static)
      for (ij = 1; ij < 13; ij++) {

        n_index = 12*index + ij;

        if (ij < 7) {
          for (i = 0; i < 3; i++) {
            P[n_index][i] = P[index][i]+ipvec[ij-1][i];
          }     
          P[n_index][4] = (r1)/(Q1)*P[index][4];
          P[n_index][3] = 1/(Q1)+P[index][3];
          BFS(P, n_index, swtch,index2);
        }
        
        else if (ij < 10) {
          for (i = 0; i < 3; i++) {
            P[n_index][i] = P[index][i]+opvec1[swtch][ij-7][i];
          }     
          P[n_index][4] = (r3)/(Q1)*P[index][4];
          P[n_index][3] = 1/(Q1)+P[index][3];
          BFS(P, n_index, (swtch+3)%4,index2);
        }
        
        else {
          for (i = 0; i < 3; i++) {
            P[n_index][i] = P[index][i]+opvec2[swtch][ij-10][i];
          }     
          P[n_index][4] = (r4)/(Q1)*P[index][4];
          P[n_index][3] = 1/(Q1)+P[index][3];
          BFS(P, n_index, (swtch+1)%4,index2);
        }
      }
    // }
  }

  else if (swtch == 1 || swtch == 3) {
    // #pragma omp parallel shared(vec1, P) private(ij)
    // {
    //   #pragma omp parallel for schedule(static)
      for (ij = 1; ij < 13; ij++) {
        
        n_index = 12*index + ij;

  	    if (ij < 7) {
  	      for (i = 0; i < 3; i++) {
  		      P[n_index][i] = P[index][i]+ipvec[ij-1][i];
  	      }		 
  	      P[n_index][4] = (r2)/(Q2)*P[index][4];
  	      P[n_index][3] = 1/(Q2)+P[index][3];
  	      BFS(P, n_index, swtch,index2);
  	    }

  	    else if (ij < 10) {
  	      for (i = 0; i < 3; i++) {
      		  P[n_index][i] = P[index][i]+opvec1[swtch][ij-7][i];
      		}		 
  	      P[n_index][4] = (r4)/(Q2)*P[index][4];
  	      P[n_index][3] = 1/(Q2)+P[index][3];
  	      BFS(P, n_index, (swtch+3)%4,index2);
  	    }

  	    else if (ij < 13) {
  	      for (i = 0; i < 3; i++) {
      		  P[n_index][i] = P[index][i]+opvec2[swtch][ij-10][i];
      		}		 
  	      P[n_index][4] = (r3)/(Q2)*P[index][4];
  	      P[n_index][3] = 1/(Q2)+P[index][3];
  	      BFS(P, n_index, (swtch+1)%4,index2);
  	    }
  	  }
    // }
  }
  return;
} 
 
int main(int argc, char** argv) {
  // omp_set_num_threads(1);
  
  clock_t start = clock();

  int i, j,k;

  long double (*P)[5];
  long double (*P2)[5];
  long double rotmf[3][3];
  long double mf[3][3];
  long double vec1[3];
  long double vec2[3];
  long int index2;
  P = malloc(P_SIZE * sizeof(long double[5]));
  P2 = malloc(P_SIZE * sizeof(long double[5]));

  int *checker;
  checker = malloc(P_SIZE * sizeof(int));

  // can't individually place after initialization
  long double cell2_duplicate[3][3] = {
            {nnd, 0*a, 0*a},
            {-nnd/2,nnd/2*sqrt(3), 0*a},
            {0*1.0, 0*1.0, 10.086*c*nnd/3.078*a/2.57218587467527*2.51866888630220}
          };
  long double ipvec_duplicate[6][3];
  long double opvec1_duplicate[4][3][3];
  long double opvec2_duplicate[4][3][3];
  long double r11 = v*exp(-E_1_l/kb/T);
  long double r22 = v*exp(-E_4_l/kb/T);
  long double r33 = v*exp(-E_1_4/kb/T);
  long double r44 = v*exp(-E_3_4/kb/T);
  long double r55 = v*exp(-E_c/kb/T);
  long double Q11 = 6.*(r11)+3.*(r33)+3.*(r44)+4.*(r55);
  long double Q22 = 6.*(r22)+3.*(r33)+3.*(r44)+4.*(r55);

  for (i =0; i< 6; i++)
    {
      for (j = 0; j < 3 ; j++)
	{
	  vec1[j] = cell2_duplicate[0][j];
	}
      mat_pow(rotmf,rotm1,i);
      mat_vec_mul(vec2,vec1,rotmf);
      for (j=0; j<3; j++)
	{
	  ipvec_duplicate[i][j] = vec2[j];
	}
    }
  for (k=0;k<4; k++)
    {
      for (i =0; i< 3; i++)
	{
	  for (j = 0; j < 3 ; j++)
	    {
	      vec1[j] = spos_Si[k][j]-spos_Si[k+1][j];
	    }
	  mat_pow(rotmf,rotm2,i);
	  mat_mul(mf,cell2_duplicate,rotmf);
	  mat_vec_mul(vec2,vec1,mf);
	  for (j=0; j<3; j++)
	    {
	      opvec1_duplicate[k][i][j] = vec2[j];
	    }
	}
    }
  for (k=0;k<4; k++)
    {
      for (i =0; i< 3; i++)
	{
	  for (j = 0; j < 3 ; j++)
	    {
	      vec1[j] = spos_Si[k+2][j]-spos_Si[k+1][j];
	    }
	  mat_pow(rotmf,rotm2,i);
	  mat_mul(mf,cell2_duplicate,rotmf);
	  mat_vec_mul(vec2,vec1,mf);
	  for (j=0; j<3; j++)
	    {
	      opvec2_duplicate[k][i][j] = vec2[j];
	    }
	}
    }
  for (i = 0; i < 3; i++) {
    memcpy(&cell2[i], &cell2_duplicate[i], sizeof(cell2_duplicate[0]));
  }
  for (i = 0; i<6; i++) {
    memcpy(&ipvec[i], &ipvec_duplicate[i], sizeof(ipvec_duplicate[0]));
  }
  for (k = 0; k < 4; k++) {
    for (i = 0; i < 3; i++) {
      memcpy(&opvec1[k][i], &opvec1_duplicate[k][i], sizeof(opvec1_duplicate[0][0]));
      memcpy(&opvec2[k][i], &opvec2_duplicate[k][i], sizeof(opvec2_duplicate[0][0]));
    }
}
  memcpy(&r1, &r11, sizeof(r11));
  memcpy(&r2, &r22, sizeof(r11));
  memcpy(&r3, &r33, sizeof(r11));
  memcpy(&r4, &r44, sizeof(r11));
  memcpy(&r5, &r55, sizeof(r11));
  memcpy(&Q1, &Q11, sizeof(r11));
  memcpy(&Q2, &Q22, sizeof(r11));
      

  r1 = v*exp(-E_1_l/kb/T);
  r2 = v*exp(-E_4_l/kb/T);
  r3 = v*exp(-E_1_4/kb/T);
  r4 = v*exp(-E_3_4/kb/T);
  r5 = v*exp(-E_c/kb/T);

  Q1 = 6*r1+3*r3+3*r4+4*r5;
  Q2 = 6*r2+3*r3+3*r4+4*r5;

  // initializes all elements to 0
  for (i = 0; i < P_SIZE; i++) {
    for (j = 0; j < 5; j++) {
      P[i][j] = 0.0;
      P2[i][j] = 0.0;
    }
  }

  /* Intializes random number generator */
  P[0][0] = 0.0;
  P[0][1] = 0.0;
  P[0][2] = 0.0;
  P[0][3] = 0.0;
  P[0][4] = 1.0;

  printf("%Lf\n", P[0][4]);

  BFS(P, 0, 0, &index2);
  printf("%ld\n",index2);
  #pragma omp parallel shared(P2, P) private(i, j)
  {
    #pragma omp parallel for schedule(static)
    for (i = 0; i < P_SIZE; i++) {
      P2[i][0] = P[i][0];
      P2[i][1] = P[i][1];
      P2[i][2] = P[i][2];
      P2[i][3] = P[i][3];
        
      for (j = i; j < P_SIZE; j++) {
        if (((checker[j] == 0 && fabsl(P[j][0] - P[i][0])<0.0000001) &&
  	   (fabsl(P[j][1] - P[i][1])<0.0000001 && fabsl(P[j][2] - P[i][2])<0.0000001)) &&
            (fabsl(P[j][3] - P[i][3])<0.0000001)) {
    	      P2[i][4] += P[j][4];
    	      checker[j] = 1;
    	  }
  	  }
       }

    for (i = 0; i < P_SIZE; i++) {
      if (P2[i][4] > 0.0) {
    	  printf("%Lf %Lf %Lf %Lf %.11Lf\n", P2[i][0], P2[i][1], P2[i][2], P2[i][3], P2[i][4]);
    	}
    }
  }

  printf("\n%f\n", (clock() - start)/(CLOCKS_PER_SEC/1000.0));

  free(P);
  free(P2);

  free(checker);

  return 0;
}

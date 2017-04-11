#include <math.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#define N 4
#define P_SIZE ((int) pow(12, N))

/* Barriers */
long double E_1_l = 2.951646;
long double E_4_l = -3590.207247+3592.956166;
long double E_1_4 = -3589.95316882+3592.97102056;
long double E_3_4 = -3589.37845071+3592.97101757;
long double E_c = 3.067482; 

/* Boltzmann constant */
long double kb = 8.6173324e-5;

/* Vibrational frequency and temperature*/
long double v = 1.6e13;
long double T = 1300;
/* T = 1800; */

/* Nearest neighbour distance */
long double nnd = 3.08472680894400e-10;

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

// mat0 is product, mat1, mat2 are things to be multiplied
void mat_mul(long double prod[3][3], long double mat1[3][3], long double mat2[3][3]) {
  int i, j, k;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      prod[i][j] = 0.0;
      for (k = 0; k < 3; k++) {
        prod[i][j] += mat1[i][k] * mat2[k][j];
      }
    }
  }

  return; 
}

// // this version tries repeated exponentiation
// void mat_pow (long double prod[3][3], long double mat[3][3], int power) {
//   int i, j, k;

//   // set to identity matrix
//   for (i = 0; i < 3; i++) {
//     for (j = 0; j < 3; j++) {
//       if (i == j) {
//         prod[i][j] = 1.0;
//       }
//       else {
//         prod[i][j] = 0.0;
//       }
//     }
//   }

//   if (power == 0) {
//     return;  
//   }
//   // repeated exponentiation
//   else {
//     long double temp[3][3];
//     while (power != 0) {
//       if (power % 2 == 1) {
//         mat_mul(temp, prod, mat);
//         for (i = 0; i < 3; i++) {
//           memcpy(&prod[i], &temp[i], sizeof(temp[0]));
//         }
//       }
      
//       power = power >> 1;
      
//       if (power == 0) {
//         break;
//       }

//       mat_mul(temp, mat, mat);
//       for (i = 0; i < 3; i++) {
//         memcpy(&mat[i], &temp[i], sizeof(temp[0]));
//       }
//     }
//   }

//   return;
// }

// void mat_vec_mul(long double prod[3], long double vec[3], long double mat[3][3]) {
//   int i, k;
  
//   // #pragma omp parallel shared(prod, vec, mat) private (i, k)
//   // {
//   //   #pragma omp for schedule(static)
//     for (i = 0; i < 3; i++) {
//       prod[i] = 0.0;
//       for (int k = 0; k < 3; k++) {
//         prod[i] += vec[k] * mat[k][i];
//       }
//     }
//   // }
  
//   return; 
// }

// mat0 is product, mat1, mat2 are things to be multiplied
// void mat_mul(long double prod[3][3], long double mat1[3][3], long double mat2[3][3]) {
//   int i, j, k;

//   // #pragma omp parallel shared(prod, mat1, mat2) private(i, j, k)
//   // {
//   //   #pragma omp for schedule(static)
//     for (i = 0; i < 3; i++) {
//       for (j = 0; j < 3; j++) {
//         prod[i][j] = 0.0;
//         for (k = 0; k < 3; k++) {
//           prod[i][j] += mat1[i][k] * mat2[k][j];
//         }
//       }
//     }
//   // }

//   return; 
// }

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
    // square mat2 and let mat0 equal mat2**2?
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          mat0[i][j] = 0.0;
          for (int k = 0; k < 3; k++) {
              mat0[i][j] += mat1[i][k]*mat2[k][j];
            }
        }
      }
    // let mat1 = mat2**2
    for (int i = 0; i < 3; i++) {
      memcpy(&mat1[i], &mat0[i], sizeof(mat0[0]));
      }
  }

  return; 
}

void mat_vec_mul(long double prod[3], long double vec[3], long double mat[3][3]) {
  int i, k;
  
  for (i = 0; i < 3; i++) {
    prod[i] = 0.0;
    for (int k = 0; k < 3; k++) {
      prod[i] += vec[k] * mat[k][i];
    }
  }
  
  return; 
}

/* Run simulation */
void BFS(long double P[P_SIZE][5], long int index, int swtch) {

  if (log(12.0*index + 1.0)/log(12.0) >= N) {
    return;
  }

  /* transition rates */
  long double r1 = v*exp(-E_1_l/kb/T);
  long double r2 = v*exp(-E_4_l/kb/T);
  long double r3 = v*exp(-E_1_4/kb/T); 
  long double r4 = v*exp(-E_3_4/kb/T); 
  long double r5 = v*exp(-E_c/kb/T);

  long double rotmf[3][3];
  long double mf[3][3];
  long double vec1[3];
  long double vec2[3];
  long double Q1,Q2;

  /* Initialize the array to store positions (llocs), which level - hexagonal
  % or cubic we are on - (k) and charge state (chg) */ 
  Q1 = 6*r1+3*r3+3*r4+4*r5;
  Q2 = 6*r2+3*r3+3*r4+4*r5;

  int i, ij;
  // new index
  int n_index;

  if (swtch == 0 || swtch == 2) {
    for (ij = 1; ij < 13; ij++) {

      n_index = 12*index + ij;

      if (ij < 7) {
          for (i = 0; i < 3; i++) {
              vec1[i] = cell2[0][i];
            }
        mat_pow(rotmf,rotm1,ij-1);
        mat_vec_mul(vec2,vec1,rotmf);
          for (i = 0; i < 3; i++) {
              P[n_index][i] = P[index][i]+vec2[i];
            }         
        P[n_index][4] = r1/Q1*P[index][4];
        P[n_index][3] = 1/Q1+P[index][3];
        BFS(P, n_index, swtch);
        }
        
      else if (ij < 10) {
          for (i = 0; i < 3; i++) {
            vec1[i] = spos_Si[swtch][i]-spos_Si[swtch+1][i];
          }
        mat_pow(rotmf,rotm2,ij-7);
        mat_mul(mf,cell2,rotmf);
        mat_vec_mul(vec2,vec1,mf);
          for (i = 0; i < 3; i++) {
              P[n_index][i] = P[index][i]+vec2[i];
            }         
        P[n_index][4] = r3/Q1*P[index][4];
        P[n_index][3] = 1/Q1+P[index][3];
        BFS(P, n_index, (swtch+3)%4);
        }
        
      else {
          for (i = 0; i < 3; i++) {
              vec1[i] = spos_Si[swtch+2][i]-spos_Si[swtch+1][i];
            }
        mat_pow(rotmf,rotm2,ij-10);
        mat_mul(mf,cell2,rotmf);
        mat_vec_mul(vec2,vec1,mf);
          for (i = 0; i < 3; i++) {
              P[n_index][i] = P[index][i]+vec2[i];
            }         
        P[n_index][4] = r4/Q1*P[index][4];
        P[n_index][3] = 1/Q1+P[index][3];
        BFS(P, n_index, (swtch+1)%4);
        }
      }
  }

  else if (swtch == 1 || swtch == 3) {
    for (ij = 1; ij < 13; ij++) {
      
      n_index = 12*index + ij;

        if (ij < 7) {
          for (i = 0; i < 3; i++) {
              vec1[i] = cell2[0][i];
            }
          mat_pow(rotmf,rotm1,ij-1);
          mat_vec_mul(vec2,vec1,rotmf);
          for (i = 0; i < 3; i++) {
              P[n_index][i] = P[index][i]+vec2[i];
          }         
          P[n_index][4] = r2/Q2*P[index][4];
          P[n_index][3] = 1/Q2+P[index][3];
          BFS(P, n_index, swtch);
        }

        else if (ij < 10) {
          for (i = 0; i < 3; i++) {
              vec1[i] = spos_Si[swtch][i]-spos_Si[swtch+1][i];
            }
          mat_pow(rotmf,rotm2,ij-7);
          mat_mul(mf,cell2,rotmf);
          mat_vec_mul(vec2,vec1,mf);
          for (i = 0; i < 3; i++) {
              P[n_index][i] = P[index][i]+vec2[i];
            }         
          P[n_index][4] = r4/Q2*P[index][4];
          P[n_index][3] = 1/Q2+P[index][3];
          BFS(P, n_index, (swtch+3)%4);
        }

        else if (ij < 13) {
          for (i = 0; i < 3; i++) {
            vec1[i] = spos_Si[swtch+2][i]-spos_Si[swtch+1][i];
          }
          mat_pow(rotmf,rotm2,ij-10);
          mat_mul(mf,cell2,rotmf);
          mat_vec_mul(vec2,vec1,mf);
          for (i = 0; i < 3; i++) {
              P[n_index][i] = P[index][i]+vec2[i];
            }         
          P[n_index][4] = r3/Q2*P[index][4];
          P[n_index][3] = 1/Q2+P[index][3];
          BFS(P, n_index, (swtch+1)%4);
        }
      }
  }
  return;
} 
 
int main(int argc, char** argv) {
  int i, j;

  long double (*P)[5];
  long double (*P2)[5];
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

  for (i = 0; i < 3; i++) {
    memcpy(&cell2[i], &cell2_duplicate, sizeof(cell2_duplicate[0]));
  }

  // initializes all elements to 0
  for (i = 0; i < P_SIZE; i++) {
    for (j = 0; j < N; j++) {
      P[i][j] = 0.0;
      P2[i][j] = 0.0;
    }
    checker[i] = 0;
  }

  /* Intializes random number generator */
  P[0][0] = 0.0;
  P[0][1] = 0.0;
  P[0][2] = 0.0;
  P[0][3] = 0.0;
  P[0][4] = 1.0;

  printf("%Lf\n", P[0][4]);

  BFS(P, 0, 0);

  for (i = 0; i < pow(12, N); i++) {
    P2[i][0] = P[i][0];
    P2[i][1] = P[i][1];
    P2[i][2] = P[i][2];
    P2[i][3] = P[i][3];
      
    for (j = i; j < pow(12, N); j++) {
        if (((checker[j] == 0 && P[j][0] == P[i][0]) &&
          (P[j][1] == P[i][1] && P[j][2] == P[i][2])) &&
          (P[j][3] == P[i][3])) {
            P2[i][4] += P[j][4];
            checker[j] = 1;
        }
      }
  }

  for (i = 0; i < pow(12, N); i++) {
    if (P2[i][4] > 0.0) {
        printf("%Lf %Lf %Lf %Lf %.11Lf\n", P2[i][0]*1e9, P2[i][1]*1e9, P2[i][2]*1e9, P2[i][3], P2[i][4]);
      }
  }

  // remember to free P, P2, checker

  return 0;
}

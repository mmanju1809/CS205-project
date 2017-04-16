#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#define N 3
#define num 5

// Energy barriers
long double E_1_l = 2.951646;
long double E_4_l = -3590.207247+3592.956166;
long double E_1_4 = -3589.95316882+3592.97102056;
long double E_3_4 = -3589.37845071+3592.97101757;
long double E_c = 3.067482; 

// Boltzmann constant
long double kb = 8.6173324e-5;

// Vibrational frequency and temperature
long double v = 1.6e13;
long double T = 1300;
/* T = 1800; */

// Nearest neighbor distance
long double nnd =  3.08472680894400e-1;

// Silicon atom positions in unit cell
long double si_cell_pos[6][3] = {
                              {0.33333333,  0.66666667,  0.93750000-1},
                              {0.00000, 0.00000, 0.1875},
                              {0.66666667,  0.33333333,  0.43750000},
                              {0.000000,  0.000000,  0.68750000},
                              {0.33333333,  0.66666667,  0.93750000},
                              {0.00000, 0.00000, 0.1875+1}
                            };

// Lattice scaling parameters and unit cell
long double a = 1.0104;
long double c = 1.0074;

// Define two rotation matrices to capture all in-plane
long double rotm_ip[3][3] = {
                            {0.5, -0.866025403784439, 0.0},
                            {0.866025403784439, 0.5, 0.0},
                            {0.0, 0.0, 1.0}
                          };

// and out-of-plane transitions
long double rotm_op[3][3] = {
                            {-0.5, -0.866025403784439, 0.0},
                            {0.866025403784439, -0.5, 0.0},
                            {0.0, 0.0, 1.0}
                          };

long double vec_ip[6][3];
long double vec_op_down[4][3][3];
long double vec_op_up[4][3][3];

void mat_mul(long double product[3][3], long double A[3][3], long double B[3][3]) {
    int i, j, k;

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            product[i][j] = 0.0;
            for (k = 0; k < 3; k++) {
                product[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return;
}

void mat_pow(long double product[3][3], long double A[3][3], int power) {
    int i, j, p;
    long double t[3][3];
    
    for (i = 0; i < 3; i++) {
        memcpy(&t[i], &A[i], sizeof(A[0]));
    }

    // if power == 0, set mat0 to identity matrix
    if (power == 0) {
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                if (i == j) {
                    product[i][j] = 1.0;
                }
                else {
                    product[i][j] = 0.0;
                }
            }
        }
    }

    // if power == 1, set mat0 to mat2
    else if (power == 1) {
        for (i = 0; i < 3; i++) {
            memcpy(&product[i], &A[i], sizeof(A[0]));
        }
    }

    for (p = 1; p < power; p++) {
        // square mat2 and let mat0 equal mat2**2
        mat_mul(product, t, A);
        // let mat1 = mat2**2
        for (i = 0; i < 3; i++) {
            memcpy(&t[i], &product[i], sizeof(product[0]));
        }
    }

    return;
}

void vec_mat_mul(long double product[3], long double vec[3], long double mat[3][3]) {
    int i, k;

    for (i = 0; i < 3; i++) {
        product[i] = 0.0;
        for (k = 0; k < 3; k++) {
            product[i] += vec[k] * mat[k][i];
        }
    }

    return;
}

// at_{k} : allowable transitions in parameter k
// p_{k} : particle value in parameter k
/*
parameters:
x : x position
y : y position
z : z position
t : time step
p : probability
v : visited
*/
void BFS(long double at_x[num][num],
          long double at_y[num][num],
          long double at_z[num][num],
          long double at_t[num][num],
          long double at_p[num][num],
          long double p_x[num],
          long double p_y[num],
          long double p_z[num],
          long double p_t[num],
          long double p_p[num],
          long double p_v[num],
          long int *counter) {
    long int lower, upper, i, j;

    // probably can cut down this loop
    for (lower = 1; lower < num; lower += 12) {
        // lower index of particles we are considering
        i = (lower - 1) / 12;
        for (j = lower; j < lower + 12; j++) {
            if (p_v[i] == 1.0) {
                p_p[j] = 1.0;
                *counter++;
                if (j % 12 == 0) {
                    p_p[i] += 1.0;
                }
                p_x[j] = at_x[i][j] + p_x[i];
                p_y[j] = at_y[i][j] + p_x[i];
                p_z[j] = at_z[i][j] + p_x[i];
                p_t[j] = at_t[i][j] + p_x[i];
                p_p[j] = at_p[i][j] * p_p[i];
            }
        }
    }
}

int main(int argc, char** argv) {
    long double r1 = v*exp(-E_1_l/kb/T);
    long double r2 = v*exp(-E_4_l/kb/T);
    long double r3 = v*exp(-E_1_4/kb/T);
    long double r4 = v*exp(-E_3_4/kb/T);
    long double r5 = v*exp(-E_c/kb/T);
    long double Q1, Q2;
    long int i, j, k, counter;

    long double *p_x, *p_y, *p_z, *p_t, *p_p, *p_v;
    long double (*at_x)[num], (*at_y)[num], (*at_z)[num], (*at_t)[num], (*at_p)[num];

    // cumulative probabilities of transitioning to particle per particle
    long double *p_cum_p;

    long double rotmf[3][3];
    long double mf[3][3];
    long double vec1[3];
    long double vec2[3];
    long int n_index;

    // transition rates
    // starting from cubic site
    Q1 = 6.*r1+3.*r3+3.*r4+4.*r5;
    // starting from hexagonal site
    Q2 = 6.*r2+3.*r3+3.*r4+4.*r5;

    p_x = malloc(num * sizeof(long double));
    p_y = malloc(num * sizeof(long double));
    p_z = malloc(num * sizeof(long double));
    p_t = malloc(num * sizeof(long double));
    p_p = malloc(num * sizeof(long double));
    p_v = malloc(num * sizeof(long double));
    at_x = malloc(num * sizeof(long double));
    at_y = malloc(num * sizeof(long double));
    at_z = malloc(num * sizeof(long double));
    at_t = malloc(num * sizeof(long double));
    at_p = malloc(num * sizeof(long double));
    p_cum_p = malloc(num * sizeof(long double));

    int *checker, *switcher;

    checker = malloc(num * sizeof(int));
    switcher = malloc(num * sizeof(int));

    // lattice vector
    long double vec_lattice[3][3] = {
        {nnd, 0*a, 0*a},
        {-nnd/2,nnd/2*sqrt(3), 0*a},
        {0*1.0, 0*1.0, 10.086*c*nnd/3.078*a/2.57218587467527*2.51866888630220}
    };

    memcpy(&vec1, &vec_lattice[0], sizeof(vec_lattice[0]));

    // in plane vector
    for (i = 0; i < 6; i++) {
        mat_pow(rotmf, rotm_ip, i);
        vec_mat_mul(vec2, vec1, rotmf);
        memcpy(&vec_ip[i], &vec2, sizeof(vec2));
    }

    // out of plane vector down
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 3; j++) {
            vec1[j] = si_cell_pos[k][j] - si_cell_pos[k+1][j];
        }
        for (k = 0; k < 3; k++) {
            mat_pow(rotmf, rotm_op, i);
            mat_mul(mf, vec_lattice, rotmf);
            vec_mat_mul(vec2, vec1, mf);
            memcpy(&vec_op_down[i][k], &vec2, sizeof(vec2));
        }
    }

    // out of plane vector up
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 3; j++) {
            vec1[j] = si_cell_pos[k][j] - si_cell_pos[k+1][j];
        }
        for (k = 0; k < 3; k++) {
            mat_pow(rotmf, rotm_op, i);
            mat_mul(mf, vec_lattice, rotmf);
            vec_mat_mul(vec2, vec1, mf);
            memcpy(&vec_op_up[i][k], &vec2, sizeof(vec2));
        }
    }

    for (i = 0; i < num; i++) {
        checker[i] = 0;
        switcher[i] = 0;
    }

    for (i = 0; i < num; i++) {
        if (switcher[i] == 0 || switcher[i] == 2) {
            for (j = 0; j < 12; j++) {
                n_index = 12*i + j + 1;
                if (n_index >= num) {
                    i = num;
                    break;
                }
                // in plane
                if (j < 6) {
                    at_x[i][n_index] = vec_ip[j][0];
                    at_y[i][n_index] = vec_ip[j][1];
                    at_z[i][n_index] = vec_ip[j][2];
                    at_t[i][n_index] = 1.0/Q1;
                    at_p[i][n_index] = r1/Q1;
                }
                // below plane
                else if (j < 9) {
                    at_x[i][n_index] = vec_op_down[switcher[i]][j-6][0];
                    at_y[i][n_index] = vec_op_down[switcher[i]][j-6][1];
                    at_z[i][n_index] = vec_op_down[switcher[i]][j-6][2];
                    at_t[i][n_index] = 1.0/Q1;
                    at_p[i][n_index] = r3/Q1;
                    switcher[n_index] = (switcher[i] + 3) % 4;
                }
                // above plane
                else {
                    at_x[i][n_index] = vec_op_up[switcher[i]][k-9][0];
                    at_y[i][n_index] = vec_op_up[switcher[i]][k-9][1];
                    at_z[i][n_index] = vec_op_up[switcher[i]][k-9][2];
                    at_t[i][n_index] = 1.0/Q1;
                    at_p[i][n_index] = r4/Q1;
                    switcher[n_index] = (switcher[i] + 1) % 4;
                }
            }
        }
        else {
            for (j = 0; j < 12; j++) {
                n_index = 12*i + j + 1;
                if (n_index >= num) {
                    i = num;
                    break;
                }
                // in plane
                if (j < 6) {
                    at_x[i][n_index] = vec_ip[j][0];
                    at_y[i][n_index] = vec_ip[j][1];
                    at_z[i][n_index] = vec_ip[j][2];
                    at_t[i][n_index] = 1.0/Q2;
                    at_p[i][n_index] = r2/Q2;
                }
                // below plane
                else if (j < 9) {
                    at_x[i][n_index] = vec_op_down[switcher[i]][j-6][0];
                    at_y[i][n_index] = vec_op_down[switcher[i]][j-6][1];
                    at_z[i][n_index] = vec_op_down[switcher[i]][j-6][2];
                    at_t[i][n_index] = 1.0/Q2;
                    at_p[i][n_index] = r4/Q2;
                    switcher[n_index] = (switcher[i] + 3) % 4;
                }
                // above plane
                else {
                    at_x[i][n_index] = vec_op_up[switcher[i]][k-9][0];
                    at_y[i][n_index] = vec_op_up[switcher[i]][k-9][1];
                    at_z[i][n_index] = vec_op_up[switcher[i]][k-9][2];
                    at_t[i][n_index] = 1.0/Q2;
                    at_p[i][n_index] = r3/Q2;
                    switcher[n_index] = (switcher[i] + 1) % 4;
                }
            }
        }
    }

    // initialize parameters for first particle
    p_x[0] = 0.0;
    p_y[0] = 0.0;
    p_z[0] = 0.0;
    p_t[0] = 0.0;
    p_p[0] = 1.0;
    p_v[0] = 1.0;
    counter = 0;

    for (i = 0; i < N; i++) {
        BFS(at_x, at_y, at_z, at_t, at_p,
            p_x, p_y, p_z, p_t, p_p, p_v,
            &counter);
    }

    for (i = 0; i < num; i++) {
        p_cum_p[i] = 0.0;

        for (j = i; j < num; j++) {
            if (((checker[j] == 0 && fabsl(p_x[j] - p_x[i]) < 0.0000001) &&
                 (fabsl(p_y[j] - p_y[i]) < 0.0000001 && fabsl(p_z[j] - p_z[i]) < 0.0000001)) &&
                 (fabsl(p_t[j] - p_t[i]) < 0.0000001)) {
                p_cum_p[i] += p_p[j];
                checker[j] = 1;
            }
        }
    }

    for (i = 0; i < num; i++) {
        if (p_cum_p[i] > 0.0) {
            printf("%Lf %Lf %Lf %Lf %.11Lf\n", p_x[i], p_y[i], p_z[i], p_t[i], p_cum_p[i]);
        }
    }

    free(p_x);
    free(p_y);
    free(p_z);
    free(p_t);
    free(p_p);
    free(p_v);
    free(at_x);
    free(at_y);
    free(at_z);
    free(at_t);
    free(at_p);
    free(p_cum_p);
    free(checker);
    free(switcher);

    return 0;
}
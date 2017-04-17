#include <math.h>

// Barriers
#define ZERO_1 = -3592.97102056;
#define E_1_l = 2.951646;
#define E_1_l_2 = -(-3584.167250 - -3581.721798);
#define E_1_l_0 = -3598.056400 - -3601.470647;
#define E_4_l = -3590.207247 - -3592.956166;
#define E_4_l_2 = -3581.852746 - -3584.077936;
#define E_4_l_0 = -3598.331221 - -3601.522984;
#define E_1_4 = -3589.95316882 - -3592.97102056;
#define E_3_4 = -3589.37845071 - ZERO_1;
#define E_3 = -3592.97101757 - ZERO_1;
#define E_3_4 = E_3_4 - E_3;
#define E_c = -3589.86261278 - -3592.97102056; // 3.100865 
#define E_1_lxp2 = -3574.65353 - -3577.708882;
#define E_1_lxm2 = -3582.279110 - -3585.078091;
#define E_1_lxm2_2 = -3573.307794 - -3575.558803;
#define E_1_lxp2_2 = -3566.993155 - -3569.576452;
#define E_1_lxp2_0 = -3582.077389 - -3585.555594;
#define E_1_lxm2_0 = -3590.988981 - -3594.273917;

// Boltzmann constant
#define kb = 8.6173324e-5;
// D = (0.0276+0.0262)/2;
// v = 109019950963228/2/pi;

// Vibrational frequency and temperature
#define v = 1.6e13;
#define T = 1300;
// T = 1800;

// Nearest neighbour distance
#define nnd = 3.08472680894400e-10;

// transition rates
float r1 = v*exp(-E_1_l/kb/T);
float r2 = v*exp(-E_4_l/kb/T);
float r3 = v*exp(-E_1_4/kb/T); // 4.8
float r4 = v*exp(-E_3_4/kb/T); // -((E_1_l+E_4_l)/2-(4.8-3.7))/k/T);%3.7
float r5 = v*exp(-E_c/kb/T);

// Silicon atom positions in unit cell
#define spos_Si = [0.33333333,  0.66666667,  0.93750000-1; 0.00000, 0.00000, 0.1875; 0.66666667,  0.33333333,  0.43750000; 0.000000,  0.000000,  0.68750000;  0.33333333,  0.66666667,  0.93750000; 0.00000, 0.00000, 0.1875+1];

// lattice scaling parameters and unit cell
#define a = 1.0104;
#define c = 1.0074;
#define cell2 = [nnd, 0*a, 0*a; -nnd/2,nnd/2*sqrt(3), 0*a; 0*1.0, 0*1.0, 10.086*c*nnd/3.078*a/2.57218587467527*2.51866888630220];

// Define two rotation matrices to capture all possible in-plane and out of plane transitions
float eul1[3] = {M_PI/3.0, 0.0, 0.0};
float eul2[3] = {2*M_PI/3.0, 0.0, 0.0};
// missing euler angle to rot matrix fxn in C

// Initialize time, number of time steps (N) and number of particles (it), 
// as well as the array to store positions (llocs), which level - hexagonal
// or cubic we are on - (k) and charge state (chg)

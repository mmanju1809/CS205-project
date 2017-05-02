#include <math.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#define N 1000
#define IT 123
#define ITER 1000

/* Barriers */
long double E_1_l = 2.951646;
long double E_4_l = -3590.207247+3592.956166;
long double E_1_4 = -3589.95316882+3592.97102056;
long double E_3_4 = -3589.37845071+3592.97101757;
long double E_c = 2.781326; //3.067482; 
long double E_cb = 2.781326 - 0.225101;

/* Boltzmann constant */
long double kb = 8.6173324e-5;

/* Vibrational frequency and temperature*/
long double v = 1.6e13;
/* T = 1800; */

/* Nearest neighbour distance */
long double nnd = 3.08472680894400e-10;

/* Silicon atom positions in unit cell */
long double spos_Si[6][3] = {{0.33333333,  0.66666667,  0.93750000-1}, {0.00000, 0.00000, 0.1875}, {0.66666667,  0.33333333,  0.43750000}, {0.000000,  0.000000,  0.68750000}, {0.33333333,  0.66666667,  0.93750000}, {0.00000, 0.00000, 0.1875+1}};

/* lattice scaling parameters and unit cell */
long double a = 1.0104;
long double c = 1.0074;

/* Define two rotation matrix to capture all possible in-plane and out of
   plane transitions */
long double rotm1[3][3] = {{0.5, -0.866025403784439, 0.0},{0.866025403784439,0.5, 0.0},{0.0,0.0,1.0}};
long double rotm2[3][3] = {{-0.5, -0.866025403784439, 0.0},{0.866025403784439,-0.5, 0.0},{0.0,0.0,1.0}};


/* Calculate potential */
long double get_pot(long double x[3],long double locs[N][3][IT],long double chg[IT],int it, int h, int y)
{ 
  /* define elementary charge (e), dielectric constant (er), and permittivity
     of free space (e0) */
  long double e = 1.60217662e-19;
  long double er = sqrt(10.03*9.66);
  long double e0 = 8.854187817e-12;
  long double eps = 2.2204e-16;
  long double pi = 3.14159265358979;
  long double inv[IT];
  long double pot = 0.0;
  /* compute distances between particles and then the potential */
  for(int i = 0; i<IT; i++)
    {
      if (i != h)
	{
	  inv[i] = sqrt(pow((x[0] - locs[y][0][i]),2)+pow((x[2] - locs[y][1][i]),2)+pow((x[3] - locs[y][2][i]),2))+eps;
	}
    }
  for(int i = 0; i<IT-1; i++)
    {
      if (i != h)
	{
	  pot += e/4/pi/e0/er*chg[i]/inv[i];
	}
    }
  return pot;
}
void mat_pow(long double mat0[3][3], long double mat2[3][3],int power)
{
  long double mat1[3][3];
  for (int i = 0; i<3;i++)
    {
      for (int j = 0; j<3; j++)
	{
	  mat1[i][j] = mat2[i][j];
	}
    }
  if (power == 0)
    {
  for (int i = 0; i< 3; i++)
    {
    for (int j = 0; j<3; j++)
      {
	if (i ==j)
	  {
	mat0[i][j] = 1.0;
	  }
	else
	  {
	mat0[i][j] = 0.0;
	  }
      }
    }
    }
  else if (power == 1)
    {
      for (int i = 0; i< 3; i++)
	{
	  for (int j = 0; j<3; j++)
	    {
	      mat0[i][j] = mat2[i][j];
	    }
	}
    }
  for (int p = 1; p< power; p++)
    {
      for (int i = 0; i <3; i++)
	{
	  for (int j = 0; j < 3; j++)
	    {
	      mat0[i][j] = 0.0;
	      for (int k = 0; k < 3; k++)
		{
		  mat0[i][j] += mat1[i][k]*mat2[k][j];
		}
	    }
	}
      for (int i = 0; i< 3; i++)
	{
	  for (int j = 0; j<3; j++)
	    {
	      mat1[i][j] = mat0[i][j];
	    }
	}
    }

  return; 
}
 
void mat_mul(long double mat0[3][3],long double mat1[3][3], long double mat2[3][3])
{
  for (int i = 0; i< 3; i++)
    {
      for (int j = 0; j<3; j++)
	{
	  mat0[i][j] = 0.0;
	  for (int k = 0; k < 3; k++)
	    {
	      mat0[i][j] += mat1[i][k]*mat2[k][j];
	    }
	}
    }
  return; 
} 

void mat_vec_mul(long double vec0[3],long double vec1[3], long double mat0[3][3])
{
  for (int i = 0; i< 3; i++)
    {
      vec0[i] = 0.0;
  for (int k = 0; k < 3; k++)
    {
      vec0[i] += vec1[k]*mat0[k][i];
    }
    }
  return; 
}

/* Run simulation */
int kMC(long double cell2[3][3], long double spos_Si[6][3], long double T,long double llocs[N][3][IT])
{ 
  /* transition rates */
  long double r1 = v*exp(-E_1_l/kb/T);
  long double r2 = v*exp(-E_4_l/kb/T);
  long double r3 = v*exp(-E_1_4/kb/T); 
  long double r4 = v*exp(-E_3_4/kb/T); 
  long double r5 = v*exp(-E_c/kb/T);
  long double r6 = v*exp(-E_cb/kb/T);
  long double rotmf[3][3];
  long double mf[3][3];
  long double vec1[3];
  long double vec2[3];
  long double locs[3];
  long double locs0[3];
  long double R[16] = {0};
  long double Q; 
  long double u1;
  long double test;
  int L;
  int E;
  int EE;
  int m;
  int h, ij,i,j,y;
  /* Initialize the array to store positions (llocs), which level - hexagonal
  % or cubic we are on - (k) and charge state (chg) */ 
  double tt = 0;
  double t; 
  int itt = IT;
  int switcher = 0;
  long double chg[IT] = {-1.0};
  int k[IT] = {0};
  int k2[IT] = {0};
  /*  llocs[0][0][0] = cell2[0][0];
      llocs[0][0][1] = -cell2[0][0];*/
  for(y = 1; y< N; y++)
    {
    t = 0;
#pragma omp parallel shared(llocs,k,k2,chg) private(h)
    {
#pragma omp parallel for schedule(static)
    for (h = 0; h< IT; h ++)
      {
	if (k[h] == 0 || k[h] == 2)
	  {
	    for (ij = 1; ij< 7; ij++)
	      {
		for (i = 0; i<3; i++)
		  {
		    vec1[i] = cell2[0][i];
		  }
		mat_pow(rotmf,rotm1,ij-1);
		mat_vec_mul(vec2,vec1,rotmf);
		for (i = 0; i<3; i++)
		  {
		    locs[i] = llocs[y-1][i][h]+vec2[i];
		    locs0[i] = llocs[y-1][i][h];
		  }		 
		R[ij] = R[ij-1]+r1*exp(-(get_pot(locs,llocs,chg,IT,h,y)-get_pot(locs0,llocs,chg,IT,h,y))/kb/T);
	      }
	    for (ij = 7; ij< 10; ij++)
	      {
		for (i = 0; i<3; i++)
		  {
		    vec1[i] = spos_Si[k[h]][i]-spos_Si[k[h]+1][i];
		  }
		mat_pow(rotmf,rotm2,ij-7);
		mat_mul(mf,cell2,rotmf);
		mat_vec_mul(vec2,vec1,mf);
		for (i = 0; i<3; i++)
		  {
		    locs[i] = llocs[y-1][i][h]+vec2[i];
		    locs0[i] = llocs[y-1][i][h];
		  }
		R[ij] = R[ij-1]+r3*exp(-(get_pot(locs,llocs,chg,IT,h,y)-get_pot(locs0,llocs,chg,IT,h,y))/kb/T);
	      }
	    for (ij = 10; ij< 13; ij++)
	      {
		for (i = 0; i<3; i++)
		  {
		    vec1[i] = spos_Si[k[h]+2][i]-spos_Si[k[h]+1][i];
		  }
		mat_pow(rotmf,rotm2,ij-10);
		mat_mul(mf,cell2,rotmf);
		mat_vec_mul(vec2,vec1,mf);
		for (i = 0; i<3; i++)
		  {
		    locs[i] = llocs[y-1][i][h]+vec2[i];
		    locs0[i] = llocs[y-1][i][h];
		  }
		R[ij] = R[ij-1]+r4*exp(-(get_pot(locs,llocs,chg,IT,h,y)-get_pot(locs0,llocs,chg,IT,h,y))/kb/T);
	      }
	    for (ij = 13; ij < 16; ij++)
	      {
		R[ij] = r5+R[ij-1];
	      }
	  }
	else if (k[h] == 1 || k[h] == 3)
	  {
	    for (ij = 1; ij< 7; ij++)
	      {
		for (i = 0; i<3; i++)
		  {
		    vec1[i] = cell2[0][i];
		  }
		mat_pow(rotmf,rotm1,ij-1);
		mat_vec_mul(vec2,vec1,rotmf);
		for (i = 0; i<3; i++)
		  {
		    locs[i] = llocs[y-1][i][h]+vec2[i];
		    locs0[i] = llocs[y-1][i][h];
		  }
		R[ij] = R[ij-1]+r2*exp(-(get_pot(locs,llocs,chg,IT,h,y)-get_pot(locs0,llocs,chg,IT,h,y))/kb/T);
	      }
	    for (ij = 7; ij< 10; ij++)
	      {
		for (i = 0; i<3; i++)
		  {
		    vec1[i] = spos_Si[k[h]][i]-spos_Si[k[h]+1][i];
		  }
		mat_pow(rotmf,rotm2,ij-7);
		mat_mul(mf,cell2,rotmf);
		mat_vec_mul(vec2,vec1,mf);
		for (i = 0; i<3; i++)
		  {
		    locs[i] = llocs[y-1][i][h]+vec2[i];
		    locs0[i] = llocs[y-1][i][h];
		  }
		R[ij] = R[ij-1]+r4*exp(-(get_pot(locs,llocs,chg,IT,h,y)-get_pot(locs0,llocs,chg,IT,h,y))/kb/T);
	      }
	    for (ij = 10; ij< 13; ij++)
	      {
		for (i = 0; i<3; i++)
		  {
		    vec1[i] = spos_Si[k[h]+2][i]-spos_Si[k[h]+1][i];
		  }
		mat_pow(rotmf,rotm2,ij-10);
		mat_mul(mf,cell2,rotmf);
		mat_vec_mul(vec2,vec1,mf);
		for (i = 0; i<3; i++)
		  {
		    locs[i] = llocs[y-1][i][h]+vec2[i];
		    locs0[i] = llocs[y-1][i][h];
		  }
		R[ij] = R[ij-1]+r3*exp(-(get_pot(locs,llocs,chg,IT,h,y)-get_pot(locs0,llocs,chg,IT,h,y))/kb/T);
	      }
	    for (ij = 13; ij < 16; ij++)
	      {
		R[ij] = r5+R[ij-1];
	      }
	  }
	else if (k[h] == 4)
	  {
            R[1] = r6;
	  }
	if (k[h] == 4)
	  {
	    Q = R[1];
	    E = 0;
	    EE = 0;
	  }
	else
	  {
	    Q = R[15];
	    E = 14;
	    EE = 14;
	  }
	u1 = (1.0*((rand()%(RAND_MAX-1))+1))/(1.0*RAND_MAX);
	test = Q*u1;
	L = 0;
	for (j = 0; j< 16; j++)
	  {
            if (L>E)
	      {
		printf("\nunsuccessful\n");
	      }
            m = (int) floor((E+L)/2);
            if (R[m]<test && R[m+1]<test)
	      {
	    L = m+1;
	      }
            else if (R[m]>test && R[m+1] > test)
	      {
		 E = m-1;
	      }
            else
	      {
	      break;
	      }
	  }
        if (m < 6 && EE == 14)
	  {
	    for (i = 0; i<3; i++)
		  {
		    vec1[i] = cell2[0][i];
		  }
	    mat_pow(rotmf,rotm1,m);
	    mat_vec_mul(vec2,vec1,rotmf);
	    for (i = 0; i<3; i++)
	      {
		llocs[y][i][h] = llocs[y-1][i][h]+vec2[i];
	      }
	  }
        else if (m <9 && EE == 14)
	  {
	    for (i = 0; i<3; i++)
	      {
		vec1[i] = spos_Si[k[h]][i]-spos_Si[k[h]+1][i];
	      }
	    mat_pow(rotmf,rotm2,m-6);
	    mat_mul(mf,cell2,rotmf);
	    mat_vec_mul(vec2,vec1,mf);
	    for (i = 0; i<3; i++)
	      {
		llocs[y][i][h] = llocs[y-1][i][h]+vec2[i];
	      }
	    k[h] = (k[h]+3)%4;
	    k2[h] = (k[h]+3)%4;
	  }
	else if (m < 12 && EE == 14)
	  {
	    for (i = 0; i<3; i++)
	      {
		vec1[i] = spos_Si[k[h]+2][i]-spos_Si[k[h]+1][i];
	      }
	    mat_pow(rotmf,rotm2,m-9);
	    mat_mul(mf,cell2,rotmf);
	    mat_vec_mul(vec2,vec1,mf);
	    for (i = 0; i<3; i++)
	      {
		llocs[y][i][h] = llocs[y-1][i][h]+vec2[i];
	      }
	    k[h] = (k[h]+1)%4;
	    k2[h] = (k[h]+1)%4;
	  }
	else if (m > 11)
	  {
	    //	    for (int i = y; i<N; i++)
	    //{
	    //	llocs[i][0][h] = llocs[y-1][0][h];
	    //	llocs[i][1][h] = llocs[y-1][1][h];
	    //	llocs[i][2][h] = llocs[y-1][2][h];
	    //}
	    k[h] = 4;
	    chg[h] = 2;
	    //if( (itt-1) < 1)
	    // { 
	    //	switcher = 1; 
	    //}  
	    itt = (itt-1) > 1? (itt-1): 1;
	  }
	else
	  {
	    k[h] = k2[h];
	    chg[h] = -1;
	    itt = (itt) >= IT? IT: (itt+1);
	    //printf("\nunsuccessful\n");
	  }
	u1 = (1.0*((rand()%(RAND_MAX-1))+1))/(1.0*RAND_MAX);
	t = t +1/Q*log(1/u1)/itt;
      }
    }
    /*    if(switcher == 1)
      {
	break;
	}*/
    #pragma omp parallel for default(shared) private(h) schedule(static) reduction(+:tt)
    for (h =0; h<IT; h++)
      {
	tt += t;
      }
    if (tt>10)
      {
	L = y;
	break;
      }

    }
  E = 0;
  for (i = 0; i < IT; i++)
    {
      if (k[i] != 4)
	{
	  printf("%Lf %Lf %Lf\n",llocs[y][0][i]*1e9,llocs[y][1][i]*1e9,llocs[y][2][i]*1e9);
	  E++;
	}
      else
	{
	  continue; 
	}
    }
  return E;


}
 
int main(int argc, char** argv)
{
  time_t t;
  long double frac; 
  int i;
  long double cell2[3][3] = {{nnd, 0*a, 0*a}, {-nnd/2,nnd/2*sqrt(3), 0*a}, {0*1.0, 0*1.0, 10.086*c*nnd/3.078*a/2.57218587467527*2.51866888630220}};
  long double T = 1023;
  /* Intializes random number generator */
  long double llocs[N][3][IT] = {0};
  long double hlocs[IT][3] = {{47.1710e-10, -38.3955e-10, -149.5400e-10}, {45.6250e-10, -37.5025e-10, -147.0100e-10}, {48.7180e-10, -35.7165e-10, -139.4100e-10}, {50.2650e-10, -38.3955e-10, -129.2900e-10}, {50.2650e-10, -34.8235e-10, -131.8200e-10}, {47.1720e-10, -36.6095e-10, -124.2200e-10}, {48.7180e-10, -37.5025e-10, -126.7500e-10}, {45.6250e-10, -32.1455e-10, -121.6900e-10}, {59.5450e-10, -20.5375e-10, -124.2200e-10}, {51.8110e-10, -14.2865e-10, -119.1600e-10}, {6.9600e-10, -13.3935e-10, -121.6900e-10}, {53.3580e-10, -9.8225e-10, -124.2200e-10}, {39.4380e-10, -30.3595e-10, -109.0400e-10}, {37.8920e-10, -27.6805e-10, -109.0400e-10}, {6.9600e-10, -20.5375e-10, -114.1000e-10}, {6.9600e-10, -18.7515e-10, -111.5700e-10}, {40.9850e-10, -16.9655e-10, -109.0400e-10}, {42.5320e-10, -10.7145e-10, -111.5700e-10}, {47.1710e-10, -6.2505e-10, -109.0400e-10}, {50.2650e-10, -8.0365e-10, -111.5700e-10}, {54.9050e-10, -7.1435e-10, -114.1000e-10}, {59.5440e-10, -8.0365e-10, -111.5700e-10}, {28.6120e-10, -24.1085e-10, -106.5000e-10}, {33.2520e-10, -23.2155e-10, -103.9700e-10}, {34.7990e-10, -24.1085e-10, -106.5000e-10}, {33.2520e-10, -21.4305e-10, -101.4400e-10}, {6.9600e-10, -13.3935e-10, -101.4400e-10}, {16.2390e-10, -13.3935e-10, -106.5000e-10}, {19.3330e-10, -13.3935e-10, -106.5000e-10}, {28.6120e-10, -13.3935e-10, -101.4400e-10}, {40.9850e-10, -6.2505e-10, -98.9100e-10}, {5.4130e-10, -62.5045e-10, -88.7800e-10}, {30.1590e-10, -30.3595e-10, -88.7800e-10}, {27.0660e-10, -26.7875e-10, -91.3200e-10}, {27.0650e-10, -19.6445e-10, -88.7800e-10}, {28.6120e-10, -18.7515e-10, -91.3200e-10}, {14.6930e-10, -16.0725e-10, -96.3800e-10}, {33.2520e-10, -12.5005e-10, -93.8500e-10}, {37.8920e-10, -8.0365e-10, -96.3800e-10}, {44.0780e-10, -2.6785e-10, -96.3800e-10}, {-42.5320e-10, -66.9695e-10, -81.1900e-10}, {11.5990e-10, -64.2905e-10, -86.2500e-10}, {-33.2520e-10, -59.8255e-10, -78.6600e-10}, {-31.7060e-10, -58.9325e-10, -81.1900e-10}, {-28.6120e-10, -53.5755e-10, -86.2500e-10}, {-27.0660e-10, -52.6825e-10, -83.7200e-10}, {-13.1460e-10, -48.2175e-10, -81.1900e-10}, {28.6120e-10, -20.5375e-10, -83.7200e-10}, {133.7820e-10, -20.5375e-10, -83.7200e-10}, {25.5190e-10, -18.7515e-10, -86.2500e-10}, {0.7730e-10, -29.4665e-10, -71.0600e-10}, {17.7860e-10, -16.0725e-10, -60.9400e-10}, {19.3330e-10, -9.8225e-10, -43.2200e-10}, {17.7860e-10, -10.7145e-10, -40.6900e-10}, {20.8790e-10, -10.7145e-10, -40.6900e-10}, {20.8790e-10, -1.7855e-10, -43.2200e-10}, {17.7860e-10, -7.1435e-10, -33.0900e-10}, {10.0530e-10, -2.6785e-10, -10.3100e-10}, {10.0530e-10, -0.8925e-10, -7.7800e-10}, {11.5990e-10, 0.0005e-10, 15.0000e-10}, {-0.7730e-10, -28.5735e-10, 78.2900e-10}, {-0.7730e-10, -26.7875e-10, 80.8200e-10}, {11.5990e-10, 0.0005e-10, 80.8200e-10}, {-22.4260e-10, -14.2865e-10, 93.4800e-10}, {-11.6000e-10, -20.5375e-10, 98.5400e-10}, {-19.3330e-10, -17.8585e-10, 98.5400e-10}, {-20.8790e-10, -18.7515e-10, 101.0700e-10}, {-17.7860e-10, -15.1795e-10, 98.5400e-10}, {-13.1460e-10, -16.0725e-10, 101.0700e-10}, {5.4130e-10, -8.9295e-10, 103.6000e-10}, {6.9590e-10, -6.2505e-10, 103.6000e-10}, {-14.6930e-10, -15.1795e-10, 108.6600e-10}, {-11.6000e-10, -13.3935e-10, 106.1300e-10}, {-10.0530e-10, -8.9295e-10, 113.7300e-10}, {-8.5070e-10, -8.0365e-10, 111.2000e-10}, {-6.9600e-10, -5.3575e-10, 111.2000e-10}, {-5.4130e-10, -2.6785e-10, 111.2000e-10}, {-2.3200e-10, -0.8925e-10, 113.7300e-10}, {-0.7730e-10, 0.0005e-10, 106.1300e-10}, {-14.6930e-10, -18.7515e-10, 121.3200e-10}, {-14.6930e-10, -15.1795e-10, 118.7900e-10}, {268.3370e-10, -292.8790e-10, 174.4800e-10}, {6.9590e-10, -6.2505e-10, 174.4800e-10}, {3.8660e-10, -0.8925e-10, 174.4800e-10}, {6.9590e-10, -0.8925e-10, 174.4800e-10}, {2.3200e-10, 0.0005e-10, 166.8900e-10}, {44.0780e-10, 0.8925e-10, -93.8500e-10}, {-70.3710e-10, 29.4665e-10, -10.3100e-10}, {-25.5190e-10, 8.9295e-10, 7.4100e-10}, {-23.9720e-10, 11.6085e-10, 7.4100e-10}, {-19.3330e-10, 14.2865e-10, 7.4100e-10}, {-20.8790e-10, 13.3945e-10, 9.9400e-10}, {-62.6380e-10, 44.6465e-10, 12.4700e-10}, {-3.8670e-10, 42.8605e-10, 4.8800e-10}, {-39.4390e-10, 47.3255e-10, 12.4700e-10}, {-40.9860e-10, 50.0035e-10, 12.4700e-10}, {10.0530e-10, 6.2505e-10, 17.5300e-10}, {8.5060e-10, 5.3575e-10, 20.0700e-10}, {8.5060e-10, 8.9295e-10, 17.5300e-10}, {17.7860e-10, 10.7155e-10, 15.0000e-10}, {-16.2400e-10, 16.0725e-10, 20.0700e-10}, {-13.1460e-10, 37.5025e-10, 20.0700e-10}, {-20.8790e-10, 43.7535e-10, 17.5300e-10}, {-34.7990e-10, 46.4325e-10, 17.5300e-10}, {-30.1590e-10, 49.1105e-10, 17.5300e-10}, {-28.6120e-10, 48.2185e-10, 15.0000e-10}, {-19.3330e-10, 48.2185e-10, 15.0000e-10}, {-31.7060e-10, 51.7895e-10, 17.5300e-10}, {-40.9850e-10, 64.2905e-10, 20.0700e-10}, {-3.8670e-10, 26.7875e-10, 30.1900e-10}, {2.3200e-10, 10.7155e-10, 40.3200e-10}, {20.8790e-10, 12.5015e-10, 42.8500e-10}, {-3.8670e-10, 12.5015e-10,52.9700e-10}, {-11.6000e-10, 40.1815e-10, 65.6300e-10}, {-19.3330e-10, 42.8605e-10, 65.6300e-10}, {-23.9730e-10, 66.9695e-10, 65.6300e-10}, {0.7730e-10, 4.4645e-10, 83.3500e-10}, {13.1460e-10, 2.6785e-10, 80.8200e-10}, {-2.3200e-10, 0.8925e-10, 139.0400e-10}, {-0.7740e-10, 1.7865e-10, 154.2300e-10}, {0.7730e-10, 0.8925e-10, 169.4200e-10}, {3.8670e-10, 0.8925e-10, 169.4200e-10}, {2.3200e-10, 1.7865e-10, 174.4800e-10}};
  for (int ii = 0; ii < IT; ii++)
    {
      llocs[0][0][ii] = hlocs[ii][0];
      llocs[0][1][ii] = hlocs[ii][1];
      llocs[0][2][ii] = hlocs[ii][2];
    }
  srand((unsigned) time(&t));
  frac = 0.0;
  #pragma omp parallel for default(shared) private(i) schedule(static) reduction(+:frac)
  for (i = 0; i<ITER; i++) 
    {
      frac += 1.0*kMC(cell2,spos_Si,T,llocs);
    }
  frac = frac/(1.0*ITER);
  printf("%Lf\n",frac);
  return 0;
}

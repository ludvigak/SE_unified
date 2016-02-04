#ifndef __SE_HASIMOTO_DECOMP
#define __SE_HASIMOTO_DECOMP

#include <math.h>

double self_coeff(double xi)
{ 
    return -4*xi/sqrt(PI); 
}

void op_A(double A[3][3], double x[3], double xi)
{
    double r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    double r = sqrt(r2);
    double b = xi*xi*r2;
    double c1 = 2*(xi*exp(-b)/(sqrt(PI)*r2) + erfc(xi*r)/(2*r*r2) );
    double c2 = r2*c1 - 4*xi*exp(-b)/sqrt(PI);
    int i, j;
    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	    A[i][j] = c1*x[i]*x[j] + c2*(i==j);
}

void op_B(double B[3][3], double k[3], double xi)
{
    double k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
    double c = k2/(4*xi*xi);
    double d = 8*PI*(1+c)*exp(-c)/(k2*k2);
    int i, j;
    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	    B[i][j] = d*( k2*(i==j) - k[i]*k[j] );
}

void op_BB(double B[3][3], double k[3], double xi)
{
    double k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
    double c = k2/(4*xi*xi);
    double d = 8*PI*(1+c)/(k2*k2);
    int i, j;
    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	    B[i][j] = d*( k2*(i==j) - k[i]*k[j] );
}

#endif

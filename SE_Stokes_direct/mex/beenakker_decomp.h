#ifndef __SE_BEENAKKER_DECOMP
#define __SE_BEENAKKER_DECOMP

#include <math.h>

double self_coeff(double xi)
{ 
    return -8*xi/sqrt(PI);
}

double C1(double r)
{
    double r2 = r*r;
    return erfc(r) + 2*(2*r2 - 3)*r*exp(-r2)/sqrt(PI);
}

double C2(double r)
{
    double r2 = r*r;
    return erfc(r) + 2*(1 - 2*r2)*r*exp(-r2)/sqrt(PI);
}

void op_A(double A[3][3], double x[3], double xi)
{
    double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    double c1 = C1(xi*r)/r;
    double c2 = C2(xi*r)/(r*r*r);
    int i,j;
    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	    A[i][j] = c1*(i==j) + c2*x[i]*x[j];
}

void op_B(double B[3][3], double k[3], double xi)
{
    double nk=sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);
    double w2 = (nk/xi)*(nk/xi);
    double k2 = nk*nk;
    double c = 8*PI*( 1.0/(w2*w2) + 1.0/(4.0*w2) + 1.0/8.0 )*
	exp(-w2/4.0)/(xi*xi*xi*xi);
    int i, j;
    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	    B[i][j] = c*( k2*(i==j) - k[i]*k[j] );
}

void op_BB(double B[3][3], double k[3], double xi)
{
    double nk=sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);
    double w2 = (nk/xi)*(nk/xi);
    double k2 = nk*nk;
    double c = 8*PI*( 1.0/(w2*w2) + 1.0/(4.0*w2) + 1.0/8.0 )/(xi*xi*xi*xi);
    int i, j;
    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	    B[i][j] = c*( k2*(i==j) - k[i]*k[j] );
}

#endif

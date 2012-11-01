#define OP_TAG "BEENAKKER FD"

// Probably excessive unrolling here, but it was fun!
// Quick benching says that 75% of op_A is spent on erfc

void op_A(double A[3][3], double x[3], double n[3], double xi)
{
    double r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    double r = sqrt(r2);
    double c = xi*xi*r2;
    double xdotn = x[0]*n[0] + x[1]*n[1] + x[2]*n[2];

    double ec = exp(-c);
    double C = -2.0/(r2*r2)*( 3.0/r*erfc(xi*r) + 2.0*xi/sqrt(PI)*(3.0+2.0*c-4.0*c*c)*ec );
    double D = 8.0/sqrt(PI)*xi*xi*xi*(2.0-c)*ec;

    A[0][0]=C*x[0]*x[0]*xdotn + D*( x[0]*n[0] + x[0]*n[0] + xdotn);
    A[0][1]=C*x[0]*x[1]*xdotn + D*( x[1]*n[0] + x[0]*n[1] );
    A[0][2]=C*x[0]*x[2]*xdotn + D*( x[2]*n[0] + x[0]*n[2] );
    A[1][1]=C*x[1]*x[1]*xdotn + D*( x[1]*n[1] + x[1]*n[1] + xdotn);
    A[1][2]=C*x[1]*x[2]*xdotn + D*( x[2]*n[1] + x[1]*n[2] );
    A[2][2]=C*x[2]*x[2]*xdotn + D*( x[2]*n[2] + x[2]*n[2] + xdotn);
    A[1][0]=A[0][1]; // Reuse symmetric part
    A[2][0]=A[0][2];
    A[2][1]=A[1][2];
}

// Calculate constants C and D for a list of distances
void op_A_CD(double* restrict C, double* restrict D, double* restrict r2l, int n, double xi)
{
    double r2,r,c;
    int i;
    for(i=0;i<n;i++) // This vectorizes
    {
	r2 = r2l[i];
	r = sqrt(r2);
	c = xi*xi*r2;
	C[i] = -2.0/(r2*r2)*( 3.0/r*erfc(xi*r) + 2.0*xi/sqrt(PI)*(3.0+2.0*c-4.0*c*c)*exp(-c) );
	D[i] = 8.0/sqrt(PI)*xi*xi*xi*(2.0-c)*exp(-c);
    }
}

// Calculate interactions between two points in both directions using stored C and D
void op_A_symm_CD(double A1[3][3], double A2[3][3], double x[3], double n1[3], double n2[3], 
		  double xi, double C, double D)
{
    // x1=x, x2=-x
    double x1dotn1 = x[0]*n1[0] + x[1]*n1[1] + x[2]*n1[2];
    double x2dotn2 = -( x[0]*n2[0] + x[1]*n2[1] + x[2]*n2[2] );
    
    for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
	{
	    A1[i][j] = C*x[i]*x[j]*x1dotn1 + D*( x[i]*n1[j] + x[j]*n1[i] + (i==j)*x1dotn1 );
	    A2[i][j] = C*x[i]*x[j]*x2dotn2 + D*(-x[i]*n2[j] - x[j]*n2[i] + (i==j)*x2dotn2 );
	}
}

// Calculate interactions between two points in both directions
void op_A_symm(double A1[3][3], double A2[3][3], double x[3], double n1[3], double n2[3], 
	       double xi, double r2)
{

    double r = sqrt(r2);
    double c = xi*xi*r2;
    double C = -2.0/(r2*r2)*( 3.0/r*erfc(xi*r) + 2.0*xi/sqrt(PI)*(3.0+2.0*c-4.0*c*c)*exp(-c) );
    double D = 8.0/sqrt(PI)*xi*xi*xi*(2.0-c)*exp(-c);

    // x1=x, x2=-x
    double x1dotn1 = x[0]*n1[0] + x[1]*n1[1] + x[2]*n1[2];
    double x2dotn2 = -( x[0]*n2[0] + x[1]*n2[1] + x[2]*n2[2] );
    
    for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
	{
	    A1[i][j] = C*x[i]*x[j]*x1dotn1 + D*( x[i]*n1[j] + x[j]*n1[i] + (i==j)*x1dotn1 );
	    A2[i][j] = C*x[i]*x[j]*x2dotn2 + D*(-x[i]*n2[j] - x[j]*n2[i] + (i==j)*x2dotn2 );
	}
}

void op_B(double Bi[3][3][3], double k[3], double xi)
{
    double k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
    double k2i = 1/k2;
    double c = k2/(xi*xi);
    double psi = -PI*(8.0+2*c+c*c)*k2i*exp(-c/4);

    Bi[0][0][0]=psi*(k[0]+k[0]+k[0]-2*k[0]*k[0]*k[0]*k2i);
    Bi[0][0][1]=psi*(k[1]-2*k[0]*k[0]*k[1]*k2i);
    Bi[0][0][2]=psi*(k[2]-2*k[0]*k[0]*k[2]*k2i);
    Bi[0][1][0]=psi*(k[1]-2*k[0]*k[1]*k[0]*k2i);
    Bi[0][1][1]=psi*(k[0]-2*k[0]*k[1]*k[1]*k2i);
    Bi[0][1][2]=psi*(-2*k[0]*k[1]*k[2]*k2i);
    Bi[0][2][0]=psi*(k[2]-2*k[0]*k[2]*k[0]*k2i);
    Bi[0][2][1]=psi*(-2*k[0]*k[2]*k[1]*k2i);
    Bi[0][2][2]=psi*(k[0]-2*k[0]*k[2]*k[2]*k2i);
    Bi[1][0][0]=psi*(k[1]-2*k[1]*k[0]*k[0]*k2i);
    Bi[1][0][1]=psi*(k[0]-2*k[1]*k[0]*k[1]*k2i);
    Bi[1][0][2]=psi*(-2*k[1]*k[0]*k[2]*k2i);
    Bi[1][1][0]=psi*(k[0]-2*k[1]*k[1]*k[0]*k2i);
    Bi[1][1][1]=psi*(k[1]+k[1]+k[1]-2*k[1]*k[1]*k[1]*k2i);
    Bi[1][1][2]=psi*(k[2]-2*k[1]*k[1]*k[2]*k2i);
    Bi[1][2][0]=psi*(-2*k[1]*k[2]*k[0]*k2i);
    Bi[1][2][1]=psi*(k[2]-2*k[1]*k[2]*k[1]*k2i);
    Bi[1][2][2]=psi*(k[1]-2*k[1]*k[2]*k[2]*k2i);
    Bi[2][0][0]=psi*(k[2]-2*k[2]*k[0]*k[0]*k2i);
    Bi[2][0][1]=psi*(-2*k[2]*k[0]*k[1]*k2i);
    Bi[2][0][2]=psi*(k[0]-2*k[2]*k[0]*k[2]*k2i);
    Bi[2][1][0]=psi*(-2*k[2]*k[1]*k[0]*k2i);
    Bi[2][1][1]=psi*(k[2]-2*k[2]*k[1]*k[1]*k2i);
    Bi[2][1][2]=psi*(k[1]-2*k[2]*k[1]*k[2]*k2i);
    Bi[2][2][0]=psi*(k[0]-2*k[2]*k[2]*k[0]*k2i);
    Bi[2][2][1]=psi*(k[1]-2*k[2]*k[2]*k[1]*k2i);
    Bi[2][2][2]=psi*(k[2]+k[2]+k[2]-2*k[2]*k[2]*k[2]*k2i);
}

void op_BB(double Bi[27], double k[3], double xi)
{
    double k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
    double k2i = 1/k2;
    double c = k2/(xi*xi);
    double psi = -PI*(8.0+2*c+c*c)*k2i;

    Bi[0]=psi*(k[0]+k[0]+k[0]-2*k[0]*k[0]*k[0]*k2i);
    Bi[1]=psi*(k[1]-2*k[0]*k[0]*k[1]*k2i);
    Bi[2]=psi*(k[2]-2*k[0]*k[0]*k[2]*k2i);
    Bi[3]=psi*(k[1]-2*k[0]*k[1]*k[0]*k2i);
    Bi[4]=psi*(k[0]-2*k[0]*k[1]*k[1]*k2i);
    Bi[5]=psi*(-2*k[0]*k[1]*k[2]*k2i);
    Bi[6]=psi*(k[2]-2*k[0]*k[2]*k[0]*k2i);
    Bi[7]=psi*(-2*k[0]*k[2]*k[1]*k2i);
    Bi[8]=psi*(k[0]-2*k[0]*k[2]*k[2]*k2i);
    Bi[9]=psi*(k[1]-2*k[1]*k[0]*k[0]*k2i);
    Bi[10]=psi*(k[0]-2*k[1]*k[0]*k[1]*k2i);
    Bi[11]=psi*(-2*k[1]*k[0]*k[2]*k2i);
    Bi[12]=psi*(k[0]-2*k[1]*k[1]*k[0]*k2i);
    Bi[13]=psi*(k[1]+k[1]+k[1]-2*k[1]*k[1]*k[1]*k2i);
    Bi[14]=psi*(k[2]-2*k[1]*k[1]*k[2]*k2i);
    Bi[15]=psi*(-2*k[1]*k[2]*k[0]*k2i);
    Bi[16]=psi*(k[2]-2*k[1]*k[2]*k[1]*k2i);
    Bi[17]=psi*(k[1]-2*k[1]*k[2]*k[2]*k2i);
    Bi[18]=psi*(k[2]-2*k[2]*k[0]*k[0]*k2i);
    Bi[19]=psi*(-2*k[2]*k[0]*k[1]*k2i);
    Bi[20]=psi*(k[0]-2*k[2]*k[0]*k[2]*k2i);
    Bi[21]=psi*(-2*k[2]*k[1]*k[0]*k2i);
    Bi[22]=psi*(k[2]-2*k[2]*k[1]*k[1]*k2i);
    Bi[23]=psi*(k[1]-2*k[2]*k[1]*k[2]*k2i);
    Bi[24]=psi*(k[0]-2*k[2]*k[2]*k[0]*k2i);
    Bi[25]=psi*(k[1]-2*k[2]*k[2]*k[1]*k2i);
    Bi[26]=psi*(k[2]+k[2]+k[2]-2*k[2]*k[2]*k[2]*k2i);

}


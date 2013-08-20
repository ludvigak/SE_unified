#define OP_TAG "BEENAKKER FD"

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

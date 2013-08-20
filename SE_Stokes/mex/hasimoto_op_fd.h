#define OP_TAG "HASIMOTO FD"

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

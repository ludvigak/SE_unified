#define OP_TAG "BEENAKKER FD"

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


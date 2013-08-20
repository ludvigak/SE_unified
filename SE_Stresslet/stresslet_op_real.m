function A = stresslet_op_real( x, nvec, xi)

r2 = x(1)^2+x(2)^2+x(3)^2;
r = sqrt(r2);
c = xi^2*r2;
xdotn = x(1)*nvec(1) + x(2)*nvec(2) + x(3)*nvec(3); 

C=-2/r2^2*( 3.0/r*erfc(xi*r) + 2.0*xi/sqrt(pi)*(3.0+2.0*c-4.0*c^2)*exp(-c) );
D=8/sqrt(pi)*xi^3*(2.0-c)*exp(-c);

A = zeros(3,3);

A(1,1)=C*x(1)*x(1)*xdotn + D*( x(1)*nvec(1) + x(1)*nvec(1) + xdotn);
A(2,1)=C*x(2)*x(1)*xdotn + D*( x(1)*nvec(2) + x(2)*nvec(1) );
A(3,1)=C*x(3)*x(1)*xdotn + D*( x(1)*nvec(3) + x(3)*nvec(1) );
A(1,2)=C*x(1)*x(2)*xdotn + D*( x(2)*nvec(1) + x(1)*nvec(2) );
A(2,2)=C*x(2)*x(2)*xdotn + D*( x(2)*nvec(2) + x(2)*nvec(2) + xdotn);
A(3,2)=C*x(3)*x(2)*xdotn + D*( x(2)*nvec(3) + x(3)*nvec(2) );
A(1,3)=C*x(1)*x(3)*xdotn + D*( x(3)*nvec(1) + x(1)*nvec(3) );
A(2,3)=C*x(2)*x(3)*xdotn + D*( x(3)*nvec(2) + x(2)*nvec(3) );
A(3,3)=C*x(3)*x(3)*xdotn + D*( x(3)*nvec(3) + x(3)*nvec(3) + xdotn);

end


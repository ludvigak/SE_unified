function A = stresslet_op_real( x, nvec, xi)
r = sqrt(x(1)^2+x(2)^2+x(3)^2);
rhat=x/r;
c = xi^2*r^2;
rhatdotn = sum(rhat.*nvec); 

psi=-2/r*( 3.0/r*erfc(xi*r) + 2.0*xi/sqrt(pi)*(3.0+2.0*c-4.0*c^2)*exp(-c) );
phi=8/sqrt(pi)*xi^3*r*(2.0-c)*exp(-c);

A=psi*rhatdotn*(rhat')*rhat+phi*(rhat'*nvec+nvec'*rhat+rhatdotn*eye(3)); 

end


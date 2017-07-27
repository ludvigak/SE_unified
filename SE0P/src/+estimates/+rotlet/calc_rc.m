function rc=calc_rc(t, opt, tol)
%Rotlet real space cut-off

xi = opt.xi;
Q = sum(t(:).^2);
V = prod(opt.box);

c = sqrt(8*Q./(3*V));
D = 4*xi^2*(c/tol)^4;

w = lambertw(D);
rc = sqrt(w)/(2*xi);

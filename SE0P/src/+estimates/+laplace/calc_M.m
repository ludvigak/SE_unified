function M = calc_M(f, opt, tol)

Q = sum(f(:).^2);
V = prod(opt.box);
xi = opt.xi;
c = sqrt(8*Q*xi^2 ./ (pi*V) );

D = 1/(3*xi^2)*(c/tol)^(4/3);
w = lambertw(D);
K = sqrt(3*xi^2*w);
M = ceil(K*opt.box(1)/pi);

if(mod(M,2))
    M = M +1;
end

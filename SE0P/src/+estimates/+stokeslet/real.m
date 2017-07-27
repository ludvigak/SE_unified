function est = real(f, opt)
% Stokeslet real space truncation erro
% est = real(f, opt);

Q = sum(f(:).^2);
V = prod(opt.box);
xi = opt.xi;
L = max(opt.box);
rc = opt.rc;

est = 2*sqrt(Q*rc/V)* exp(-(rc*xi)^2);

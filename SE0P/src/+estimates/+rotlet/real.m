function est = real(f, opt)
% Rotlet real space truncation error
% est = real(f, opt);

rc = opt.rc;
xi = opt.xi;
Q = sum(f(:).^2);
V = prod(opt.box);

est = sqrt(8*Q./(3*V*rc)).*exp(-(xi*rc).^2);

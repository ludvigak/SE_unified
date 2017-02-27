function est = stokeslet_est_k(f, opt, xi)

K = pi*opt.M(1)/opt.box(1);
Q = sum(f(:).^2);
V = prod(opt.box);
L = max(opt.box);
kinf = L*K/2/pi;

est = sqrt(8*Q*K/(3*pi*L^3*xi^2)) * exp(-K^2/4/xi^2);
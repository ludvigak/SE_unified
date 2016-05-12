function est = stokeslet_est_k_Lindbo(f, opt, xi)

K = pi*opt.M(1)/opt.box(1);
Q = sum(f(:).^2);
V = prod(opt.box);
L = max(opt.box);
kinf = L*K/2/pi;

est = sqrt(Q) * xi^3*L^2/pi^4./kinf.^(3/2) .* exp(-(pi*kinf/xi/L).^2);

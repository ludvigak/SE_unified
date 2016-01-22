function est = rotlet_est_k(t, opt, xi)

K = pi*opt.M(1)/opt.box(1);
Q = sum(t(:).^2);
V = prod(opt.box);

est = sqrt(8*xi^2*Q./(3*pi*V*K)) .* exp(-K.^2/4/xi^2);

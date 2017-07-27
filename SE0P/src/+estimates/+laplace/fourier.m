function est = fourier(f, opt)

K = pi*opt.M(1)/opt.box(1);
Q = sum(f(:).^2);
V = prod(opt.box);
xi = opt.xi;

est = sqrt(8*Q*xi^2 ./ (pi*V*K.^3) ) .* exp(-K.^2/4/xi^2);

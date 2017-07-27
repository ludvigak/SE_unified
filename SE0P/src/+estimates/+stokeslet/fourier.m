function est = fourier(f, opt)
% Stokeslet fourier space truncation error
% est = fourier(f,opt)

K = pi*opt.M(1)/opt.box(1);
Q = sum(f(:).^2);
V = prod(opt.box);
xi = opt.xi;
L = max(opt.box);
R = opt.R;

est = sqrt(Q) * R * K^3 / (pi * xi^2 * L) * exp(-K^2/4/xi^2);

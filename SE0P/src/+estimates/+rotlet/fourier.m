function est = fourier(f, opt)
% Rotlet fourier space truncation error
% est = fourier(f,opt)

K = pi*opt.M(1)/opt.box(1);
Q = sum(f(:).^2);
V = prod(opt.box);

est = sqrt(8*opt.xi^2*Q./(3*pi*V*K)) .* exp(-K.^2/4/opt.xi^2);

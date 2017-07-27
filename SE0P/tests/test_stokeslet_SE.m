clear
% Test stokeslet Spectral Ewald

rng(1);
L = 1;
box = [L L L];
N = 100;
[x, f] = vector_system(N, box, 3);
x(1,:) = 0;
x(2,:) = box;
M0 = 10; 
opt.M = M0*box;
opt.xi = pi*M0 / 12;
opt.P = 24;
opt.oversampling = 4;
opt.rc = 6 / opt.xi;
opt.box = box;

% Ewald
pre = stokeslet_precomp(opt);
uf = stokeslet_fourier_space(x, f, opt, pre);
ur = stokeslet_real_space(x, f, opt);
us = -4*opt.xi/sqrt(pi)*f;
ue = uf+ur+us;
% Direct
u = stokeslet_direct(x, f, box);

error = norm(u-ue, inf) / norm(u,inf);
assert(error < 1e-13, 'stokeslet SE did not match direct sum')
clear
% Test stresslet Spectral Ewald

rng(1);
L = 1;
box = [L L L];
N = 100;
[x, n] = vector_system(N, box, 3);
[~, q] = vector_system(N, box, 3);

x(1,:) = 0;
x(2,:) = box;
M0 = 10; 
opt.M = M0*box;
opt.xi = pi*M0 / 14;
opt.P = 32;
opt.oversampling = 1+sqrt(3);
opt.rc = 6 / opt.xi;
opt.box = box;

% Ewald
f = [n q];
pre = stresslet_precomp(opt);
uf = stresslet_fourier_space(x, f, opt, pre);
ur = stresslet_real_space(x, f, opt);
ue = uf + ur;
% Direct
u = stresslet_direct(x, f, box);

error = norm(u-ue, inf) / norm(u,inf)
assert(error < 1e-13, 'stresslet SE did not match direct sum')
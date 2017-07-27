clear
rng(1);

L = 1;
box = [L L L];
N = 1000;
[x, f] = vector_system(N, box, 3);

% Maxmimize point-to-point distance
x(1,:) = 0;
x(2,:) = box;

M0 = 30; % Set M0, the rest is auto

opt.M = M0*box;
opt.xi = pi*M0 / 12;
opt.P = 24;
opt.oversampling = 1+sqrt(3);
opt.rc = 6 / opt.xi;
opt.box = box;
opt

% Precomp
tic
pre = stokeslet_precomp(opt);
toc
% Ewald
tic
uf = stokeslet_fourier_space(x, f, opt, pre);
ur = stokeslet_real_space(x, f, opt);
us = -4*opt.xi/sqrt(pi)*f;
ue = uf+ur+us;
toc
% > 1.3

% Direct
tic
u = stokeslet_direct(x, f, box);
toc

error = norm(u-ue, inf) / norm(u,inf)
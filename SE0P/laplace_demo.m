clear
rng(1);

N = 10000;
L = 3;
box = [L L L];
[x, f] = laplace_system(N, box);

M0 = 28; % Set M0, the rest is auto

opt.M = M0*box;
opt.xi = pi*M0 / 12;
opt.P = 24;
opt.oversampling = 3.2;
opt.rc = 6 / opt.xi;
opt.box = box

% Precomp
disp('Precomputing...')
tic
pre = laplace_precomp(opt);
toc


% Ewald
disp('Free-space Ewald...')
tic
uf = laplace_fourier_space(x, f, opt, pre);
ur = laplace_real_space(x, f, opt);
us = -f*opt.xi *2/sqrt(pi);
ue = uf+ur+us;
toc
% Direct
disp('Direct sum...')
tic
u = laplace_direct(x, f, box);
toc

error = norm(u-ue, inf) / norm(u,inf)
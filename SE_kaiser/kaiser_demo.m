clf

rng(1)
box = [1 1 1];   % domain
N = 1000;         % number of charged particles

M0 = 28;

opt.M = M0*box;
opt.xi = pi*M0 / 12;opt.xi=4;
opt.rc = 6 / opt.xi;
opt.box = box;

opt.layers = (opt.M(3)-1)/2;

opt.beta = 2*pi*.98^2;
p = [0.0059524 -0.19643 6.1905];
opt.beta = polyval(p,opt.xi);

% charge-neutral system
[x q] = SE_charged_system(N,box,'scalar');

idx = 1:10;
% compute FD Ewald sum
ref = SE3P_direct_fd_mex(idx,x,q,opt);

rms = @(x) sqrt(sum(x.^2)/numel(x));

opt.P = 16;
[u time]= spectral_ewald_kaiser(1:N,x,q,opt);

e_rms   = rms(u(idx)-ref)/rms(ref)

clc; clear all

restoredefaultpath
init

%tol = 1e-12;

rng(2)
box = [1 1 1];   % domain
N = 100;          % number of charged particles

M0 = 28;

opt.M = M0*box;
opt.xi = pi*M0 / 12;
opt.rc = 6 / opt.xi;
opt.box = box;

opt.layers = (opt.M(1)-1)/2;

opt.beta = 1.4*pi*.98^2;

% charge-neutral system
[x q] = vector_system(N,box);

idx = 1:min(N,10);
% compute FD Ewald sum
ref1p = SE1P_direct_fd_mex(idx,x,q,opt) + SE1P_direct_k0_mex(idx,x,q,opt);
ref2p = SE2P_direct_fd_mex(idx,x,q,opt) + SE2P_direct_k0_mex(idx,x,q,opt);
ref3p = SE3P_direct_fd_mex(idx,x,q,opt);

rms = @(x) sqrt(sum(x.^2)/numel(x));

opt.P = 14;
opt.M = M0*box;
[u time] = se3p_fourier_space_kaiser(1:N,x,q,opt);
kaiser_3p = rms(u(idx)-ref3p)/rms(ref3p);

opt.M = M0*box(1)-6;
opt.s = 3.5;
opt.s0= 2.5;
opt.n = 5;
[u time]= se2p_fourier_space_kaiser(1:N,x,q,opt);
kaiser_2p   = rms(u(idx)-ref2p)/rms(ref2p);

opt.M = M0*box(1);
opt.s = 3.5;
opt.s0= 2.5;
opt.n = 5;
[u time] = se1p_fourier_space_kaiser(1:N,x,q,opt);
kaiser_1p= rms(u(idx)-ref1p)/rms(ref1p);

addpath('../SE')
opt.P = 20;
opt.M = M0*box;
u    = spectral_ewald(1:N,x,q,opt.xi,opt);
gaussian_3p   = rms(u(idx)-ref3p)/rms(ref3p);

addpath('../SE2P/src')
opt.M = M0*box(1);
u    = se2p_fourier_space(x,q,opt);
gaussian_2p   = rms(u(idx)-ref2p)/rms(ref2p);

addpath('../SE1P/src')
opt.M = M0*box(1);
opt.sl=opt.s;
opt.nl=5;
u      = se1p_fourier_space(x,q,opt);
gaussian_1p   = rms(u(idx)-ref1p)/rms(ref1p);

var={'Function','Three_P','Two_P','One_P'};
table({'KAISER','GAUSSIAN'}',[kaiser_3p gaussian_3p]',[kaiser_2p gaussian_2p]',[kaiser_1p gaussian_1p]','variablenames',var)
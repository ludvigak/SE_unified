% 2P Spectral Ewald, basic accuracy/convergence computation

clear all,  close all

rng(1)
N = 20;
L = 4;
box = [L L L];
[x, f] = vector_system(N, box);
M0 = 20;		  % Set M0 to an even number, the rest is auto

opt.M = M0*box(1);
opt.xi = pi*M0 / 12;
opt.P = 32;
opt.rc = 6 / opt.xi;
opt.box = box;
opt.layers = 40;
opt.s = 4;			% upsamling rate on local pad
opt.n=8;			% local pad
opt.s0= 2;			% upsampling rate of zero mode

% % compute FD Ewald sum
ref  = SE2P_direct_fd_mex(1:N,x,f,opt);
ref0 = se2p_k0_direct(1:N,x,f,opt,false);
ref=ref+ref0;

u = se2p_fourier_space(x,f,opt);
% compute RMS error
rms_err = rms(u-ref)/rms(ref)

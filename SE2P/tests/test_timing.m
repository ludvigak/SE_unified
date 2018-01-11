% 2P Spectral Ewald, timing comparison with and without adaptive FFT and
% different upsampling rates. 

% TOL = 1e-12

clear all,  close all

rng(1)
N = 400000;
L = 8;
box = [L L L];
[x, f] = vector_system(N, box);
M0 = 20; % Set M0 to an even number, the rest is auto

opt.M = M0*box(1);
opt.xi = pi*M0 / 12;
opt.P = 24;
opt.box = box;
opt.s = 4;
opt.s0= 4;
% opt.n= NOT defined and therefore is set to the maximum possible value.

% compute reference
[ref t_full]= se2p_fourier_space(x,f,opt);

% compute solution with double precision
opt.P = 24;
opt.s=4;
opt.s0=2;
opt.n=22;

[u t_opt]= se2p_fourier_space(x,f,opt);

t_opt_upsampling = t_opt.total;
t_full_upsampling = t_full.total;

% compute RMS error
rms_err = rms(u-ref)/rms(ref);

table(t_full_upsampling,t_opt_upsampling,rms_err, 'variableNames',{'t_full_upsampling' 't_opt_upsampling','rms_err'})
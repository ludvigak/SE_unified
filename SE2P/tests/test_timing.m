% 2P Spectral Ewald, basic accuracy/convergence computation

clear all,  close all

rng(1)
N = 200000;
L = 1;
box = 8*[1 1 1];
[x, f] = vector_system(N, box);
M0 = 14; % Set M0 to an even number, the rest is auto

opt.M = M0*box(1);
opt.xi = pi*M0 / 12;
opt.P = 32;
opt.rc = 6 / opt.xi;
opt.box = box;
opt.layers = 50;
opt.s = 4;
opt.n=58;
opt.s0= 2;

% compute reference
[ref t_full]= se2p_fourier_space(x,f,opt);

% compute solution with double precision
opt.P = 24;
opt.n=48;
[u t_opt]= se2p_fourier_space(x,f,opt);

t_opt_upsampling = t_opt.total;
t_full_upsampling = t_full.total;

% compute RMS error
rms_err = rms(u-ref)/rms(ref);

table(t_full_upsampling,t_opt_upsampling,rms_err, 'variableNames',{'t_full_upsampling' 't_opt_upsampling','rms_err'})

clear
rng(1);

N = 50;
L = 1;
box = [L L L];
[x, f] = vector_system(N, box);

M0 = 20; % Set M0 to an even number, the rest is auto

opt.M = M0*box;
opt.xi = pi*M0 / 12;
opt.P = 32;
opt.box = box;
opt.layers = 10;


%% Ewald
disp('Three periodic Ewald...')
uf= se3p_fourier_space(x, f, opt);

%% Direct
t=tic;
ref = SE3P_direct_fd_mex(1:N,x,f,opt);
tdirect=toc(t);
rms_error = rms(ref-uf) / rms(ref)

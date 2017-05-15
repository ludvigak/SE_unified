clear
rng(1);

N = 20;
L = 1;
box = [L L L];
%box = [1.3 1.2 1];
[x, f] = vector_system(N, box);

M0 = 20; % Set M0 to an even number, the rest is auto

opt.M = M0*box(3);
opt.xi = pi*M0 / 12;
opt.P = 28;
opt.rc = 6 / opt.xi;opt.rc = min(opt.rc,L);
opt.box = box;
opt.layers = 100;
opt.sl = 2.5;
opt.nl = 4;
opt.s0 = 2.5;


%% Ewald
disp('One periodic Ewald...')
[uf tfourier]= se1p_fourier_space(x, f, opt);
[ur treal]= se1p_real_space(1:N,x, f, opt,false);
us = -f*opt.xi *2/sqrt(pi);
ue = uf+ur+us;

%% Direct
disp('Direct sum...')
% 1/r decays very slowly. Do not use the se1p_direct unless you
% know what is happening
%[u tdirect]= se1p_direct(1:N, x, f, box,opt);
t=tic;
u1 = SE1P_direct_rsrc_mex(1:N,x,f,opt);
u2 = SE1P_direct_fd_mex(1:N,x,f,opt);
u3 = SE1P_direct_k0_mex(1:N,x,f,opt);
u4 = SE1P_direct_self_mex(1:N,x,f,opt);
tdirect=toc(t);
u = u1+u2+u3+u4;
rms_error = rms(u-ue) / rms(u)

fprintf('fourier-sp time\t real-sp time\t total time\t direct time\n');
fprintf('%10.6f\t%10.6f\t%10.6f\t%10.6f\n',...
       tfourier.total,treal,tfourier.total+treal,tdirect);
clear all, close all,

N = 100000;
SE_opt.box=[1 1 1];
SE_opt.P = 16;
SE_opt.M = 30; 
SE_opt.sampling_factor=2;
SE_opt.zLim = [-.5 1.5];
xi = 8;

[x f] = SE_charged_system(N,SE_opt.box,'vector');

tic
u = SE2P_Stokes(1:N,x,f,xi,SE_opt);
toc

tic
SE_static = SE2P_Stokes_pre(x,xi,SE_opt);
toc

tic
[uu stats] = SE2P_Stokes(1:N,x,f,xi,SE_opt,SE_static);
toc

max(abs(u(:)-uu(:)))
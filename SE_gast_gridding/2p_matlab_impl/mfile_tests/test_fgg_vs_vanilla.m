clear all, close all,

N = 100;
rand('state',1)

box = [1 1 1];
opt.M = 31;
opt.Mz=2*opt.M;
opt.P = 19;
opt.h=box(1)/opt.M;

w=opt.h*(opt.P-1)/2;
opt.a=-w;

x = rand(N,3);
q = rand(N,1);

c=10;
opt.c=10;
opt.box=box;

% grid
tic, H  = SE2P_grid(x,q,opt,c);            toc
tic, Hfgg = SE_fg_grid_mex(x,q,opt);   toc
max(abs(H(:)-Hfgg(:)))

% int
tic, phi = SE2P_int(x,H,opt,c);            toc
tic, phi_fgg = SE_fg_int_mex(x,H,opt); toc
max(abs(phi-phi_fgg))


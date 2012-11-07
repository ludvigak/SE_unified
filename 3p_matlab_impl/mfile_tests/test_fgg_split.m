clear all

rand('state',11)

box = [1 1 1];
M = 60*[1 1 1];
P = 16;

N = 100000;
c = 100;

opt.c=c;
opt.M=M;
opt.box=box;
opt.P=P;

x = rand(N,3);
q = rand(N,1);
disp('PRE-comp')



tic
[zx zy zz idx] = SE_fgg_expand_all_mex(x,opt);

[~, s] = sort(idx);
x = x(s,:);
q = q(s);
[zx zy zz idx] = SE_fgg_expand_all_mex(x,opt);
toc

tic
zs = SE_fgg_base_gaussian_mex(opt);
toc

disp('FG GRID')

tic
H = SE_fg_grid_mex(x,q,opt);
toc

tic
HH = SE_fg_grid_split_mex(x,q,opt,zs,zx,zy,zz,idx);
toc

disp('FG INT')

tic
p = SE_fg_int_mex(x,H,opt);
toc

tic
pp = SE_fg_int_split_mex(x,H,opt,zs,zx,zy,zz,idx);
toc

disp('ERR')
max(abs(H(:)-HH(:)))
max(abs(p-pp))

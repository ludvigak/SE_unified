function SE_static  = SE_FGG_precomp(x,xi,opt)
% SPECTRAL EWALD 3P, pre-computation of FGG vectors
% Fast Ewald method for electrostatic potential calculation, k-space part.

verb = false;

% get parameters
opt=parse_params(opt);
eta=(2*opt.w*xi/opt.m)^2;
opt.c = 2*xi^2/eta;
opt.box=opt.box;
opt.a=opt.wbox(3,1);

[zx zy zz idx] = SE_fgg_expand_all_mex(x,opt);

[~, s] = sort(idx);
x = x(s,:);

[zx zy zz idx] = SE_fgg_expand_all_mex(x,opt);
zs = SE_fgg_base_gaussian_mex(opt);

SE_static.zs=zs;
SE_static.zx=zx;
SE_static.zy=zy;
SE_static.zz=zz;
SE_static.idx=idx;
SE_static.perm=s';
SE_static.iperm(s)=1:length(s);
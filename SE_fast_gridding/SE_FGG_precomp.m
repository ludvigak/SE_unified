function [SE_static walltime]  = SE_FGG_precomp(x,xi,opt)
% SPECTRAL EWALD 3P, pre-computation of FGG vectors
% Fast Ewald method for electrostatic potential calculation, k-space part.

verb = false;

% get parameters

x = recenter_points(x, opt.box);

[zx zy zz idx] = SE_fgg_expand_all_mex(x,opt);

[idx, s] = sort(idx);
x = x(s,:);
zx = zx(:,s);
zy = zy(:,s);
zz = zz(:,s);
t=tic;
zs = SE_fgg_base_gaussian_mex(opt);
walltime=toc(t);

SE_static.zs=zs;
SE_static.zx=zx;
SE_static.zy=zy;
SE_static.zz=zz;
SE_static.idx=idx;
SE_static.perm=s';
SE_static.iperm(s)=1:length(s);
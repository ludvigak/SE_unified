function [SE_static pre] = se0p_precomp_kaiser(x,opt)
% SPECTRAL EWALD 0P, pre-computation of FKG vectors
% Fast Ewald method for electrostatic potential calculation, k-space part.

verb = false;
x = recenter_points(x, opt.box);

[zx zy zz idx] = SE_fkg_expand_all_mex(x,opt);
%[zx zy zz idx] = SE_fgg_expand_all_mex(x,opt);

[idx, s] = sort(idx);
x = x(s,:);
zx = zx(:,s);
zy = zy(:,s);
zz = zz(:,s);

%zs = SE_fgg_base_gaussian_mex(opt);
%SE_static.zs=zs;

SE_static.zx=zx;
SE_static.zy=zy;
SE_static.zz=zz;
SE_static.idx=idx;
SE_static.perm=s';
SE_static.iperm(s)=1:length(s);
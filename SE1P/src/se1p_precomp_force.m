function SE_static  = se1p_precomp_force(x,opt)
% SPECTRAL EWALD 1P, pre-computation of FGG vectors
% Fast Ewald method for electrostatic potential calculation, k-space part.

verb = false;

% get parameters
%opt=se1p_parse_params(opt);

[zx, zy, zz, zfx, zfy, zfz, idx] = SE_fgg_expand_all_force_mex_1p(x,opt);

[~, s] = sort(idx);
x = x(s,:);

[zx, zy, zz, zfx, zfy, zfz, idx] = SE_fgg_expand_all_force_mex_1p(x,opt);
zs = SE_fgg_base_gaussian_mex_1p(opt);

SE_static.zs=zs;
SE_static.zx=zx;
SE_static.zy=zy;
SE_static.zz=zz;
SE_static.zfx=zfx;
SE_static.zfy=zfy;
SE_static.zfz=zfz;
SE_static.idx=idx;
SE_static.perm=s';
SE_static.iperm(s)=1:length(s);

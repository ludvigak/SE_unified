function SE_static  = se2p_precomp(x,opt) 
% SPECTRAL EWALD 2P, pre-computation of FGG vectors
% Fast Ewald method for electrostatic potential calculation, k-space part.

verb = false;

[zx zy zz idx] = SE_fgg_expand_all_mex_2p(x,opt);

[~, s] = sort(idx);
x = x(s,:);

[zx zy zz idx] = SE_fgg_expand_all_mex_2p(x,opt);
zs = SE_fgg_base_gaussian_mex_2p(opt);

SE_static.zs=zs;
SE_static.zx=zx;
SE_static.zy=zy;
SE_static.zz=zz;
SE_static.idx=idx;
SE_static.perm=s';
SE_static.iperm(s)=1:length(s);
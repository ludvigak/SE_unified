function SE_static  = se3p_precomp_force(x,opt) 

[zx, zy, zz, zfx, zfy, zfz, idx] = SE_fgg_expand_all_force_mex(x,opt);

[~, s] = sort(idx);
x = x(s,:);

[zx, zy, zz, zfx, zfy, zfz, idx] = SE_fgg_expand_all_force_mex(x,opt);
zs = SE_fgg_base_gaussian_mex(opt);

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
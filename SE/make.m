cc = 'gcc';

cflags = '-std=c99 -fPIC -msse2';

switch cc
    case 'icc'
        coptimflags = '-O3 -static -xHOST -vec-report';
        ldoptimflags = '-O3 -static -xHOST -vec-report';
    case 'gcc'
        coptimflags = '-Wall -O3 -ffast-math';
        ldoptimflags = '-O3';
end
cdebugflags='';
lddebugflags='';

CC = ['CC=''' cc ''''];
CFLAGS = [' CFLAGS=''' cflags ''''];
DEBUGFLAGS = [' CDEBUGFLAGS=''' cdebugflags '''' ' LDDEBUGFLAGS=''' lddebugflags ''''];
OPTIMFLAGS = [' COPTIMFLAGS=''' coptimflags '''' ' LDOPTIMFLAGS=''' ldoptimflags ''''];

mex_string = ['mex ' CC CFLAGS DEBUGFLAGS OPTIMFLAGS ' -outdir bin/'];

% build FGG code from ../SE_fast_gridding
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fg_grid_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fg_int_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fgg_expand_all_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fgg_base_gaussian_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fg_grid_split_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fg_int_split_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])

% build direct summation code from ../SE_direct
eval([mex_string ' -DVERBOSE -DEWALD_REAL -DTHREE_PERIODIC ../SE_direct/SE_direct_mex.c ../SE_direct/SE3P_direct.c -o SE3P_direct_real_mex'])
eval([mex_string ' -DVERBOSE -DEWALD_RSRC -DTHREE_PERIODIC ../SE_direct/SE_direct_mex.c ../SE_direct/SE3P_direct.c -o SE3P_direct_rsrc_mex'])
eval([mex_string ' -DVERBOSE -DEWALD_FD -DTHREE_PERIODIC ../SE_direct/SE_direct_mex.c ../SE_direct/SE3P_direct.c -o SE3P_direct_fd_mex'])
eval([mex_string ' -DVERBOSE -DEWALD_SELF -DTHREE_PERIODIC ../SE_direct/SE_direct_mex.c ../SE_direct/SE3P_direct.c -o SE3P_direct_self_mex'])

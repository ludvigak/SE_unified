% cc = 'gcc';
cc = 'icc';

cflags = '-std=c99 -fPIC -msse2';

switch cc
    case 'icc'
%         coptimflags = '-O3 -static -xHOST -vec-report';
%         ldoptimflags = '-O3 -static -xHOST -vec-report';

%          coptimflags = '-O3 -static -xSSE4.1 -vec-report';
%          ldoptimflags = '-O3 -static -xSSE4.1 -vec-report';

       coptimflags = '-O3 -static -vec-report'; % -opt-report
       ldoptimflags = '-O3 -static -vec-report';
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

eval([mex_string ' -largeArrayDims -DVERBOSE -DBEENAKKER mex/stresslet_real_rc.c'])
return

% build FGG code from ../SE_fast_gridding
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fg_grid_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fg_int_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fgg_expand_all_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fgg_base_gaussian_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fg_grid_split_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fg_int_split_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])

% build direct summation code from ../SE_Stokes_direct
eval([mex_string ' -DVERBOSE -DEWALD_REAL -DTHREE_PERIODIC -DHASIMOTO ../SE_Stokes_direct/SE_Stokes_direct_mex.c ../SE_Stokes_direct/SE3P_Stokes_direct.c -o SE3P_Stokes_direct_real_mex'])
eval([mex_string ' -DVERBOSE -DEWALD_RSRC -DTHREE_PERIODIC -DHASIMOTO ../SE_Stokes_direct/SE_Stokes_direct_mex.c ../SE_Stokes_direct/SE3P_Stokes_direct.c -o SE3P_Stokes_direct_rsrc_mex'])
eval([mex_string ' -DVERBOSE -DEWALD_FD -DTHREE_PERIODIC -DHASIMOTO ../SE_Stokes_direct/SE_Stokes_direct_mex.c ../SE_Stokes_direct/SE3P_Stokes_direct.c -o SE3P_Stokes_direct_fd_mex'])
eval([mex_string ' -DVERBOSE -DEWALD_SELF -DTHREE_PERIODIC -DHASIMOTO ../SE_Stokes_direct/SE_Stokes_direct_mex.c ../SE_Stokes_direct/SE3P_Stokes_direct.c -o SE3P_Stokes_direct_self_mex'])

% build SE_Stresslet specific mex stuff
eval([mex_string ' -DVERBOSE -DBEENAKKER mex/stresslet_fast_k_scaling.c'])
eval([mex_string ' -DBEENAKKER mex/stresslet_direct_real_mexcore.c'])
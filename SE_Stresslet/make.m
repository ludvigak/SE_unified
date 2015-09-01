cc = 'gcc';
%cc = 'icc';

openmp = true;


cflags = '-std=c99 -fPIC -mavx2 -mfma';
ldflags = ' -lm ';

switch cc
    case 'icc'
%        coptimflags = '-O3 -static -vec-report -xHOST'; % -opt-report
%        ldoptimflags = '-O3 -static -vec-report -xHOST';

       coptimflags = '-O3 -static -vec-report'; % -opt-report
       ldoptimflags = '-O3 -static -vec-report';

       if openmp
            ldflags = [ldflags ' -openmp '];
            cflags = [cflags ' -openmp '];
       end
    case 'gcc'
        coptimflags = '-Wall -O3 -ffast-math -ftree-vectorizer-verbose=0';        
        ldoptimflags = '-O3';
        if openmp
            ldflags = [ldflags ' -fopenmp '];
            cflags = [cflags ' -fopenmp '];
        end
        cc = 'gcc-5';
end
cdebugflags='';
lddebugflags='';

CC = ['CC="' cc '"'];
CFLAGS = [' CFLAGS="' cflags '"'];
LDFLAGS = [' LDFLAGS="\$LDFLAGS ' ldflags  '" '];
DEBUGFLAGS = [' CDEBUGFLAGS='''  cdebugflags '''' ' LDDEBUGFLAGS=''' lddebugflags ''''];
OPTIMFLAGS = [' COPTIMFLAGS=''' coptimflags '''' ' LDOPTIMFLAGS=''' ldoptimflags ''''];

mex_string = ['mex ' CC CFLAGS LDFLAGS DEBUGFLAGS OPTIMFLAGS ' -outdir bin/'];

% build FGG code from ../SE_fast_gridding
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fg_grid_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fg_int_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fgg_expand_all_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fgg_base_gaussian_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fg_grid_split_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fg_int_split_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])

% Not doing Stokes here anymore
% % build direct summation code from ../SE_Stokes_direct
% eval([mex_string ' -DVERBOSE -DEWALD_REAL -DTHREE_PERIODIC -DHASIMOTO ../SE_Stokes_direct/SE_Stokes_direct_mex.c ../SE_Stokes_direct/SE3P_Stokes_direct.c -output SE3P_Stokes_direct_real_mex'])
% eval([mex_string ' -DVERBOSE -DEWALD_RSRC -DTHREE_PERIODIC -DHASIMOTO ../SE_Stokes_direct/SE_Stokes_direct_mex.c ../SE_Stokes_direct/SE3P_Stokes_direct.c -output SE3P_Stokes_direct_rsrc_mex'])
% eval([mex_string ' -DVERBOSE -DEWALD_FD -DTHREE_PERIODIC -DHASIMOTO ../SE_Stokes_direct/SE_Stokes_direct_mex.c ../SE_Stokes_direct/SE3P_Stokes_direct.c -output SE3P_Stokes_direct_fd_mex'])
% eval([mex_string ' -DVERBOSE -DEWALD_SELF -DTHREE_PERIODIC -DHASIMOTO ../SE_Stokes_direct/SE_Stokes_direct_mex.c ../SE_Stokes_direct/SE3P_Stokes_direct.c -output SE3P_Stokes_direct_self_mex'])

% build SE_Stresslet specific mex stuff
eval([mex_string ' -DVERBOSE -DBEENAKKER mex/stresslet_fast_k_scaling.c'])
eval([mex_string ' -DBEENAKKER mex/stresslet_direct_real_mexcore.c'])
eval([mex_string ' -largeArrayDims -DVERBOSE -DBEENAKKER mex/stresslet_real_rc.c mex/cell_list.c mex/stresslet_real_rc_mex.c -output stresslet_real_rc_mex'])
eval([mex_string ' -largeArrayDims -DVERBOSE -DBEENAKKER mex/stresslet_real_rc.c mex/cell_list.c mex/stresslet_real_rc_nomatrix_mex.c -output stresslet_real_rc_nomatrix_mex'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE mex/stresslet_fg_grid_mex.c -I../SE_fast_gridding/ ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])

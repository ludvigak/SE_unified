cc = 'icc';

cflags = '-std=c99 -fPIC ';

switch cc
    case 'icc'
        coptimflags = '-O3 -xHOST -vec-report';
        ldoptimflags = '-O3 -xHOST -vec-report';
        ompflag = ' -openmp ';
    case 'gcc'
        coptimflags = '-O3';
        ldoptimflags = '-O3';
        ompflag = ' -fopenmp ';
end
cdebugflags='';
lddebugflags='';

CC = ['CC=''' cc ''''];
CFLAGS = [' CFLAGS=''' cflags ''''];
DEBUGFLAGS = [' CDEBUGFLAGS=''' cdebugflags '''' ' LDDEBUGFLAGS=''' lddebugflags ''''];
OPTIMFLAGS = [' COPTIMFLAGS=''' coptimflags '''' ' LDOPTIMFLAGS=''' ldoptimflags ''''];

mex_string = ['mex ' CC CFLAGS DEBUGFLAGS OPTIMFLAGS];

eval([mex_string ' -DVERBOSE -DEWALD_REAL -DTHREE_PERIODIC -DHASIMOTO ' ...
    '../SE_Stokes_direct_mex.c ../SE3P_Stokes_direct.c -o SE3P_Stokes_direct_real_mex'])
eval([mex_string ' -DVERBOSE -DEWALD_RSRC -DTHREE_PERIODIC -DHASIMOTO ' ...
    '../SE_Stokes_direct_mex.c ../SE3P_Stokes_direct.c -o SE3P_Stokes_direct_rsrc_mex'])
eval([mex_string ' -DVERBOSE -DEWALD_FD -DTHREE_PERIODIC -DHASIMOTO ' ...
    '../SE_Stokes_direct_mex.c ../SE3P_Stokes_direct.c -o SE3P_Stokes_direct_fd_mex'])
eval([mex_string ' -DVERBOSE -DEWALD_SELF -DTHREE_PERIODIC -DHASIMOTO ' ...
      '../SE_Stokes_direct_mex.c ../SE3P_Stokes_direct.c -o SE3P_Stokes_direct_self_mex'])

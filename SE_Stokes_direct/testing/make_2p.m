cc = 'gcc';

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

eval([mex_string ' -DVERBOSE -DEWALD_REAL -DTWO_PERIODIC -DHASIMOTO ' ...
    '../SE_Stokes_direct_mex.c ../SE2P_Stokes_direct.c -output SE2P_Stokes_direct_real_mex'])
eval([mex_string ' -DVERBOSE -DEWALD_RSRC -DTWO_PERIODIC -DHASIMOTO ' ...
    '../SE_Stokes_direct_mex.c ../SE2P_Stokes_direct.c -output SE2P_Stokes_direct_rsrc_mex'])
eval([mex_string ' -DVERBOSE -DEWALD_FD -DTWO_PERIODIC -DHASIMOTO ' ...
    '../SE_Stokes_direct_mex.c ../SE2P_Stokes_direct.c -output SE2P_Stokes_direct_fd_mex'])
eval([mex_string ' -DVERBOSE -DEWALD_K0 -DTWO_PERIODIC -DHASIMOTO ' ...
    '../SE_Stokes_direct_mex.c ../SE2P_Stokes_direct.c -output SE2P_Stokes_direct_k0_mex'])
eval([mex_string ' -DVERBOSE -DEWALD_SELF -DTWO_PERIODIC -DHASIMOTO ' ...
      '../SE_Stokes_direct_mex.c ../SE2P_Stokes_direct.c -output SE2P_Stokes_direct_self_mex'])

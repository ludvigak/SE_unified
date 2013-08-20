cc = 'gcc';

cflags = '-std=c99 -fPIC -msse4.2';

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

mex_string = ['mex ' CC CFLAGS DEBUGFLAGS OPTIMFLAGS];

eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../mex/SE_fg_grid_mex.c ../SE_fgg.c ../SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../mex/SE_fg_int_mex.c ../SE_fgg.c ../SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../mex/SE_fgg_expand_all_mex.c ../SE_fgg.c ../SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../mex/SE_fgg_base_gaussian_mex.c ../SE_fgg.c ../SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../mex/SE_fg_grid_split_mex.c ../SE_fgg.c ../SE_fgg_MEX_params.c'])
eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../mex/SE_fg_int_split_mex.c ../SE_fgg.c ../SE_fgg_MEX_params.c'])
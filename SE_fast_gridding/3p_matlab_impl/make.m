cc = 'gcc';

cflags = '-std=c99 -fPIC -mavx2 -mfma';

switch cc
    case 'icc'
        coptimflags = '-O3 -static -xHOST -vec-report -openmp';
        ldoptimflags = '-O3 -static -xHOST -vec-report -openmp';
    case 'gcc'
        coptimflags = '-Wall -O3 -ffast-math -fopenmp';
        ldoptimflags = '-O3 -fopenmp';
        cc = 'gcc-5';
end
cdebugflags='';
lddebugflags='';

CC = ['CC=''' cc ''''];
MEXFLAGS = ' -I../ -DTHREE_PERIODIC -DVERBOSE ';
CFLAGS = [' CFLAGS=''' cflags ''''];
DEBUGFLAGS = [' CDEBUGFLAGS=''' cdebugflags '''' ' LDDEBUGFLAGS=''' lddebugflags ''''];
OPTIMFLAGS = [' COPTIMFLAGS=''' coptimflags '''' ' LDOPTIMFLAGS=''' ldoptimflags ''''];

mex_string = ['mex ' MEXFLAGS CC CFLAGS DEBUGFLAGS OPTIMFLAGS ' -outdir bin/ '];

eval([mex_string '../mex/SE_fg_grid_mex.c ../SE_fgg.c ../SE_fgg_MEX_params.c'])
eval([mex_string '../mex/SE_fg_int_mex.c ../SE_fgg.c ../SE_fgg_MEX_params.c'])
eval([mex_string '../mex/SE_fgg_expand_all_mex.c ../SE_fgg.c ../SE_fgg_MEX_params.c'])
eval([mex_string '../mex/SE_fgg_base_gaussian_mex.c ../SE_fgg.c ../SE_fgg_MEX_params.c'])
eval([mex_string '../mex/SE_fg_grid_split_mex.c ../SE_fgg.c ../SE_fgg_MEX_params.c'])
eval([mex_string '../mex/SE_fg_int_split_mex.c ../SE_fgg.c ../SE_fgg_MEX_params.c'])

eval([mex_string '-DFGG_THRD ../mex/SE_fg_grid_split_thrd_mex.c ../SE_fgg.c ../SE_fgg_MEX_params.c'])
eval([mex_string '-DFGG_THRD ../mex/SE_fg_grid_thrd_mex.c ../SE_fgg.c ../SE_fgg_MEX_params.c'])

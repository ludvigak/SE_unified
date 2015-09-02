cc = 'gcc';
%cc = 'icc';

openmp = true;


cflags = '-std=c99 -fPIC -mavx2 -mfma';
ldflags = ' -lm';

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
        coptimflags = '-Wall -O3 -ffast-math';        
        ldoptimflags = '-O3';
        if openmp
            ldflags = [ldflags ' -fopenmp '];
            cflags = [cflags ' -fopenmp '];
        end
        cc = 'gcc-5';
end
cdebugflags='';
lddebugflags='';

INC = ' -I../SE_fast_gridding  ';
CC = ['CC="' cc '"'];
CFLAGS = [' CFLAGS="' cflags '"'];
LDFLAGS = [' LDFLAGS="\$LDFLAGS ' ldflags  '" '];
DEBUGFLAGS = [' CDEBUGFLAGS='''  cdebugflags '''' ' LDDEBUGFLAGS=''' lddebugflags ''''];
OPTIMFLAGS = [' COPTIMFLAGS=''' coptimflags '''' ' LDOPTIMFLAGS=''' ldoptimflags ''''];

mex_string = ['mex ' INC CC CFLAGS LDFLAGS DEBUGFLAGS OPTIMFLAGS ' -outdir bin/'];

eval([mex_string ' -DFGG_THRD -DTHREE_PERIODIC -DVERBOSE SE_fg_grid_thrd_mex.c SE_fgg_thrd.c ' ...
      '../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c' ...
      ' ' ...
     ])

eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fg_grid_split_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])

eval([mex_string ' -DTHREE_PERIODIC -DVERBOSE ../SE_fast_gridding/mex/SE_fg_grid_mex.c ../SE_fast_gridding/SE_fgg.c ../SE_fast_gridding/SE_fgg_MEX_params.c'])

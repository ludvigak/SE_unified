% download stfmmlib at http://www.cims.nyu.edu/cmcl/fmm3dlib/fmm3dlib.html

nsource = 4000;
do_direct = true;
do_fmm    = true;
do_fse    = true;

% do sum zero
U_FMM = 0; fmm_time = 0; 
U_FSE = 0;  fse_time= 0;  
direct_time_their   = 0;
direct_time_our     = 0;

% vector system
rng(1)
L = 5.4288;
box = L*[1 1 1];
fse = fsewald('stokeslet');
fse.real_sum = @stokeslet_real_space_mkl;
[x f] = fse.generator(nsource, box);

source   = x';
ifsingle = 1;
sigma_dl = f';
ifdouble = 0;
sigma_dv = f';
sigma_sl = f';

% no need to target since source and target are the same.
target = source(:,1:nsource);
target(1,:) = target(1,:);
[ndim,ntarget] = size(target);
ntarget = 0;

% just do the potential at the source points
ifpot=1;
ifgrad=0;
ifpottarg = 0;
ifgradtarg = 0;

%% FMM
fprintf('Stokes particle target FMM in R^3\n')

%set the accuracy
iprec=1;

if(do_fmm)
    t1=tic;
    [U_FMM]=stfmm3dpart(iprec,nsource,source,ifsingle,sigma_sl, ...
                        ifdouble,sigma_dl,sigma_dv,ifpot,ifgrad,...
                        ntarget,target,ifpottarg,ifgradtarg);
    fmm_time=toc(t1);
    U_FMM = U_FMM.pot; 
end

%% FMM direct
fprintf('Stokes particle direct evaluation in R^3\n')
if(do_direct)
    t2=tic;
    [U_direct_their]=st3dpartdirect(nsource,source,ifsingle, ...
                                    sigma_sl,ifdouble,sigma_dl,...
                                    sigma_dv,ifpot,ifgrad,ntarget, ...
                                    target,ifpottarg,ifgradtarg);
    direct_time_their=toc(t2);
    %U_direct = U_direct_their;
    %dlmwrite('data/U_direct.dat',U_direct.pot, 'precision',16);
    U_direct_their = U_direct_their.pot;
else
    U_direct_their = dlmread('data/U_direct.dat');
    U_direct_their = reshape(U_direct_their,nsource,3)';
end


%% FSE
fprintf('Stokes particle FSE evaluation in R^3 ours\n')
if(do_fse)
    box = L*[1 1 1];
    xi = 6.3;
    rc = .46;
    M = 100;
    opt.M = M*[1 1 1];
    opt.xi = xi;
    opt.P = 10;
    opt.oversampling = 1+sqrt(3);
    opt.rc = rc;
    opt.box = box;
    opt.no_extra_support = true;

    % Precomp
    t = tic;
    pre = fse.precomp(opt);
    fse_pre = toc(t);
    t = tic();
    [uf wf] = fse.fourier_sum(x, f, opt, pre);
    fse_fourier = toc(t);
    fse_fft = wf.fft;
    fse_grid = wf.grid;
    fse_int = wf.int;
    fse_scale = wf.scale;
    t = tic();
    ur = fse.real_sum(x, f, opt);
    fse_real = toc(t);
    t = tic();
    us = fse.self(f, opt);
    fse_self = toc(t);

    fse_comp = fse_real + fse_fourier + fse_self;
    fse_time = fse_comp + fse_pre;
    U_FSE = (uf+ur+us);
    U_FSE = U_FSE'/2;  % difference to the direct from stfmmlib.
end

%% FSE direct
fprintf('Stokes particle direct evaluation in R^3 ours\n')
if(do_direct)
    t5 = tic;
    U_direct_our = fse.direct_sum(x,f,box);
    direct_time_our = toc(t5);
    U_direct_our = U_direct_our';
else
    U_direct_our = dlmread('data/U_direct.dat');
    U_direct_our = reshape(U_direct_our,nsource,3)';
end

table(fmm_time, fse_time, direct_time_their, direct_time_our)
table(fse_pre, fse_fft, fse_grid, fse_int, fse_scale, fse_real, fse_self)


% compare against the their direct solution, with maxtrix norm2 !!!
% as in stfmmlib
rms_error_fmm = norm((U_FMM - U_direct_their),2)/sqrt(nsource);
rel_error_fmm = norm((U_FMM - U_direct_their),2)/norm(U_direct_their, 2);

rms_error_fse = norm((U_FSE - U_direct_their),2)/sqrt(nsource);
rel_error_fse = norm((U_FSE - U_direct_their),2)/norm(U_direct_their,2);

table(rms_error_fmm, rel_error_fmm, rms_error_fse, rel_error_fse)

% compare against the their direct solution, with norm2
rms_error_fmm = norm((U_FMM(:) - U_direct_their(:)),2)/sqrt(nsource);
rel_error_fmm = norm((U_FMM(:) - U_direct_their(:)),2)/norm(U_direct_their(:), 2);

rms_error_fse = norm((U_FSE(:) - U_direct_their(:)),2)/sqrt(nsource);
rel_error_fse = norm((U_FSE(:) - U_direct_their(:)),2)/norm(U_direct_their(:),2);

table(rms_error_fmm, rel_error_fmm, rms_error_fse, rel_error_fse)

% with rms norm
rms = @(u) sqrt(1/nsource*sum(u(:).^2));
rms_error_fmm = rms(U_FMM - U_direct_their);
rel_error_fmm = rms(U_FMM - U_direct_their)/rms(U_direct_their);

rms_error_fse = rms(U_FSE - U_direct_their);
rel_error_fse = rms(U_FSE - U_direct_their)/rms(U_direct_their);

table(rms_error_fmm, rel_error_fmm, rms_error_fse, rel_error_fse)

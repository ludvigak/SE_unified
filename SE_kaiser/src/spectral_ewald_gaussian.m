function [phi walltime] = spectral_ewald_gaussian(eval_idx,x,q,opt)

fgg = true;

verb = true;
% initialize time array
walltime = struct('pre',0,'grid',0,'fft',0,'ifft',0,'scale',0,'int',0,'pregauss',0);

% parameters and constants
opt = parse_params(opt);

%gridding precomp
if(fgg)
    pre_t = tic;
    [S pt]= SE_FGG_precomp(x,opt.xi,opt);
    grid_fcn = @(f) SE_fg_grid_split_thrd_mex(x(S.perm,:),f(S.perm),opt,...
					S.zs,S.zx,S.zy,S.zz, S.idx);
    walltime.pre = toc(pre_t)-pt;
    walltime.pregauss = pt;
    % Integrator
    SI = S;
    iperm = @(u) u(SI.iperm,:);
    int_fcn = @(F) iperm(SE_fg_int_split_mex(0,F,opt,SI.zs,SI.zx,SI.zy,SI.zz,SI.idx));
else
    grid_fcn = @(f) SE_fg_grid_mex(x,f,opt);
    int_fcn = @(f) SE_fg_int_mex(x,f,opt);
end

% to grid function
grid_t=tic;
H = grid_fcn(q);
walltime.grid = toc(grid_t);

% transform and shift
fft_t=tic;
H = fftn(H);
walltime.fft = toc(fft_t);

% k-vectors
[k1 k2 k3] = k_vectors(opt.M,opt.box); 
[K1 K2 K3] = ndgrid(k1,k2,k3);

% scale
k2 = K1.^2 + K2.^2 + K3.^2;

scale_t=tic;
Z = exp(-(1-opt.eta)*k2/(4*opt.xi^2))./k2;
Z(1,1,1) = 0;
H = H.*Z;
walltime.scale = toc(scale_t);

% inverse shift and inverse transform
ifft_t=tic;
H = ifftn( H );
walltime.ifft = toc(ifft_t);

% spread and integrate
int_t=tic;
phi = int_fcn(H);
walltime.int = toc(int_t);

phi  = 4*pi*phi(eval_idx);
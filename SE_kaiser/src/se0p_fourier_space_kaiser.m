function [u varargout]= se0p_fourier_space_kaiser(x, f, opt)

% initialize time array
walltime = struct('pre',0,'prefft',0,'grid',0,'fft',0,'scale',0,'int',0);

% Check
assert(isfield(opt, 'xi'), 'xi must be given in opt struct')
assert(isfield(opt, 's'), 'oversampling rate must be given in opt struct')

% Setup vars, modify opt
fse_opt = se0p_parse_params(opt);

% precompute free-space kernel
pre_t = tic;
pre_kernel = se0p_kernel_precomp(fse_opt);
walltime.prefft = toc(pre_t);


x = bsxfun(@plus, x, fse_opt.delta);

opt = fse_opt;
opt.box = fse_opt.extended_box; % Work with box extended for Gaussian support
opt.M = fse_opt.extended_M;

% precompute kaiser window: This substitutes exp(eta*k^2/4xi^2)
pre_t = tic;
pre_window = se0p_window_precomp(fse_opt);
walltime.prefft = walltime.prefft + toc(pre_t);

% === Use vectorized code
% Gridder
% precompute zx, zy, zz / and or zs
pre_t = tic;
S = se0p_precomp_kaiser(x,opt);
walltime.pre = toc(pre_t);
grid_fcn = @(f) SE_fg_grid_split_kaiser_mex(x(S.perm,:),f(S.perm),opt,...
					   S.zx,S.zy,S.zz, S.idx);

% Integrator
SI = S;
iperm = @(u) u(SI.iperm,:);
int_fcn = @(F) iperm(SE_fg_int_split_kaiser_mex(0,F,opt,SI.zx,SI.zy,SI.zz,SI.idx));

% === Uncomment for direct code
%grid_fcn = @(f) SE_fg_grid_kaiser_mex(x,f,opt);
%int_fcn = @(f) SE_fg_int_kaiser_mex(x,f,opt);

% to grid
grid_t = tic;
H = grid_fcn(f);
walltime.grid = walltime.grid + toc(grid_t);

% transform and shift, let FFT do padding
fft_t = tic;
H = fftn(H, fse_opt.padded_M);
walltime.fft = walltime.fft + toc(fft_t);

% scale
scale_t = tic;
G = se0p_k_scaling_kaiser(H, fse_opt, pre_kernel, pre_window);
walltime.scale = toc(scale_t);


% inverse shift and inverse transform
fft_t = tic;
G = ifftn( G);
walltime.fft = walltime.fft + toc(fft_t);

% Option 1: Truncate grid before integration
M = fse_opt.extended_M;

% Option 2: Integrate directly on padded grid
%opt.M = fse_opt.padded_M;
%opt.box = fse_opt.padded_box;

int_t = tic;
u = int_fcn(G(1:M(1), 1:M(2), 1:M(3)));
walltime.int = walltime.int + toc(int_t);

if nargout==2
    walltime.total = sum(struct2array(walltime));
    varargout{1} = walltime;
end


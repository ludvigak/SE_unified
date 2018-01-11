function [u varargout]= fse_fourier_space(x, f, opt, pre, k_scaling)

% initialize time array
walltime = struct('pre',0,'grid',0,'fft',0,'scale',0,'int',0);

% Check
assert(isfield(opt, 'xi'), 'xi must be given in opt struct')
assert(isfield(opt, 'oversampling'), 'oversampling rate must be given in opt struct')

% Setup vars, modify opt
fse_opt = setup_fse(opt);
x = bsxfun(@plus, x, fse_opt.delta);
opt.box = fse_opt.extended_box; % Work with box extended for Gaussian support
opt.M = fse_opt.extended_M;

opt = parse_params(opt);
fsize = size(f);
N = fsize(1);
dim_in = fsize(2:end);

% === Use vectorized code
% Gridder
pre_t = tic;
S = SE_FGG_precomp(x,opt.xi,opt);
walltime.pre = toc(pre_t);
grid_fcn = @(f) SE_fg_grid_split_thrd_mex(x(S.perm,:),f(S.perm),opt,S.zs,S.zx,S.zy,S.zz, ...
					 S.idx);
% Integrator
SI = S;
iperm = @(u) u(SI.iperm,:);
int_fcn = @(F) iperm(SE_fg_int_split_mex(0,F,opt,SI.zs,SI.zx,SI.zy,SI.zz,SI.idx));

% === Uncomment for direct code
%grid_fcn = @(f) SE_fg_grid_mex(x,f,opt);
%int_fcn = @(f) SE_fg_int_mex(x,f,opt);


% grid + FFT
H = cell([dim_in, 1]);
for i=1:prod(dim_in)
    grid_t = tic;
    % Grid
    tmp = grid_fcn(f(:,i));
    walltime.grid = walltime.grid + toc(grid_t);
    % transform and shift, let FFT do padding
    fft_t = tic;
    H{i} = fftshift( fftn(tmp, fse_opt.padded_M) );
    walltime.fft = walltime.fft + toc(fft_t);
end

% scale
scale_t = tic;
G = k_scaling(H, fse_opt, opt, pre);
walltime.scale = toc(scale_t);
dim_out = numel(G);

% inverse shift and inverse transform
for i=1:dim_out
    fft_t = tic;
    G{i} = ifftn( ifftshift(G{i}));
    walltime.fft = walltime.fft + toc(fft_t);
end

% Option 1: Truncate grid before integration
M = fse_opt.extended_M;

% Option 2: Integrate directly on padded grid
%opt.M = fse_opt.padded_M;
%opt.box = fse_opt.padded_box;

u = zeros(N, dim_out);
for i=1:dim_out
    int_t = tic;
    u(:,i) = int_fcn(G{i}(1:M(1), 1:M(2), 1:M(3)));
    walltime.int = walltime.int + toc(int_t);
end

if nargout==2
    walltime.total = sum(struct2array(walltime));
    varargout{1} = walltime;
end


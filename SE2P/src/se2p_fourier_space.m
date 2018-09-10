function [u varargout]= se2p_fourier_space(x, f, opt)
% Fast Ewald method for electrostatic potential calculation, k-space part.

% initialize time array
walltime = struct('pre',0,'grid',0,'fft',0,'scale',0,'int',0);

% Check
assert(isfield(opt, 'xi'), 'xi must be given in opt struct')
assert(isfield(opt, 's'), 'oversampling rate must be given in opt struct')
assert(isfield(opt, 's0'), 'zero mod oversampling rate must be given in opt struct')

% parameters and constants
se_opt = se2p_parse_params(opt);


fsize = size(f);
N = fsize(1);

% === Use vectorized code
% Gridder
pre_t = tic;
S = se2p_precomp(x,se_opt);
walltime.pre = toc(pre_t);
grid_fcn = @(f) SE_fg_grid_split_thrd_mex_2p(x(S.perm,:),f(S.perm), ...
					    se_opt,S.zs,S.zx,S.zy,S.zz,S.idx);
% Integrator
SI = S;
iperm = @(u) u(SI.iperm,:);
int_fcn = @(F) iperm(SE_fg_int_split_mex_2p(0,F,se_opt,... 
					   SI.zs,SI.zx,SI.zy,SI.zz,SI.idx));

% === Uncomment for direct code
% grid_fcn = @(f) SE_fg_grid_mex_2p(x,f,se_opt);
% int_fcn = @(f) SE_fg_int_mex_2p(x,f,se_opt);

% Grid
grid_t = tic;
H = grid_fcn(f);
walltime.grid = walltime.grid + toc(grid_t);

% transform and shift
fft_t = tic;
[G Gr G0] = fftnd(H,se_opt);
walltime.fft = walltime.fft + toc(fft_t);

% scale 
scale_t = tic;
[G Gr G0] = se2p_k_scaling(G,Gr,G0,se_opt);
walltime.scale = toc(scale_t);

% inverse shift and inverse transform
fft_t = tic;
F = ifftnd(G,Gr,G0,se_opt);
walltime.fft = walltime.fft + toc(fft_t);

int_t = tic;
u = 4*pi*int_fcn(real(F));
walltime.int = walltime.int + toc(int_t);

if nargout==2
    walltime.total = sum(struct2array(walltime));
    varargout{1} = walltime;
end
function [u varargout]= se2p_fourier_space(x, f, opt)
% Fast Ewald method for electrostatic potential calculation, k-space part.

% Check
assert(isfield(opt, 'xi'), 'xi must be given in opt struct')
assert(isfield(opt, 's'), 'oversampling rate must be given in opt struct')

% parameters and constants
se_opt = se2p_parse_params(opt);


fsize = size(f);
N = fsize(1);

% === Use vectorized code
% Gridder
pre_t = tic;
S = se2p_precomp(x,se_opt);
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
H = grid_fcn(f);

% transform and shift
%G = fftn(H, [se_opt.M se_opt.M se_opt.Mz*se_opt.s]);

[G G0] = fftnd(H,se_opt.M,se_opt.M,se_opt.Mz,se_opt.s0,se_opt.s);

[G G0] = se2p_k_scaling(G,G0,se_opt);

% inverse shift and inverse transform
%F = ifftn( G ); % this should be real to eps accuracy!
%F = F(:,:,1:se_opt.Mz);
F = ifftnd(G,G0,se_opt.M,se_opt.M,se_opt.Mz,se_opt.s0,se_opt.s);

u = 4*pi*int_fcn(F);
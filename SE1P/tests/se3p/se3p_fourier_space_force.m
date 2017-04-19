function [u varargout]  = se3p_fourier_space_force(x, q, opt)

% initialize time array
walltime = struct('pre',0,'grid',0,'fft',0,'scale',0,'int',0);

assert(isfield(opt, 'xi'), 'xi must be given in opt struct')
assert(isfield(opt, 'M'), 'M must be given in opt struct')

% parameters and constants
se_opt = se3p_parse_params(opt);


%Gridder
pre_t = tic;
S = se3p_precomp_force(x,se_opt);
walltime.pre = toc(pre_t);
grid_fcn = @(f) SE_fg_grid_split_mex(x(S.perm,:),f(S.perm), ...
                                     se_opt,S.zs,S.zx,S.zy,S.zz,S.idx);
% % % Integrator
SI = S;
iperm = @(u) u(SI.iperm,:);
int_fcn = @(F) iperm(SE_fg_int_split_force_mex(x,F,se_opt,... 
                                               SI.zs,SI.zx,SI.zy,SI.zz,...
                                               SI.zfx, SI.zfy, SI.zfz, ...
                                               SI.idx));

% === Uncomment for direct code
%grid_fcn = @(F) SE_fg_grid_force_mex(x,q,se_opt);
%int_fcn = @(F) SE_fg_int_force_mexp(x,q,F,se_opt);


grid_t = tic;
%grid 
H = grid_fcn(q);
walltime.grid = walltime.grid + toc(grid_t);

% transform and shift
fft_t = tic;
G = fftn(H);
walltime.fft = walltime.fft + toc(fft_t);

% scaling
scale_t = tic;
G = se3p_k_scaling(G,se_opt);
walltime.scale = toc(scale_t);

fft_t = tic;
% inverse shift and inverse transform
F = ifftn( G );
walltime.fft = walltime.fft + toc(fft_t);

u = zeros(numel(q), 3);
int_t = tic;
u = 4*pi*int_fcn(F);
N = size(x,1);

walltime.int = walltime.int + toc(int_t);

if nargout==2
    walltime.total = sum(struct2array(walltime));
    varargout{1} = walltime;
end
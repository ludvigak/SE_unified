function [u varargout]= se1p_fourier_space_force(x, f, opt) 

% initialize time array
walltime = struct('pre',0,'grid',0,'fft',0,'scale',0,'int',0);

% Check
assert(isfield(opt, 'xi'), 'xi must be given in opt struct')
assert(isfield(opt, 'sl'), 'oversampling rate must be given in opt struct')
assert(isfield(opt, 's0'), 'k0 oversampling rate must be given in opt struct')
assert(isfield(opt, 'nl'), 'number of k0 oversamped modes must be given in opt struct')

% Setup vars, modify opt

se_opt = se1p_parse_params(opt);

fsize = size(f);
N = fsize(1);

% === Use vectorized code
% Gridder
pre_t = tic;
S = se1p_precomp_force(x,se_opt);
walltime.pre = toc(pre_t);
grid_fcn = @(f) SE_fg_grid_split_thrd_mex_1p(x(S.perm,:),f(S.perm), ...
                                             se_opt,S.zs,S.zx,S.zy,S.zz,S.idx);
% % % Integrator
SI = S;
iperm = @(u) u(SI.iperm,:);
int_fcn = @(F) iperm(SE_fg_int_split_force_mex_1p(x,F,se_opt,... 
                                                  SI.zs,SI.zx,SI.zy,SI.zz,...
                                                  SI.zfx, SI.zfy, SI.zfz, ...
                                                  SI.idx));

% === Uncomment for direct code

%grid_fcn = @(F) SE_fg_grid_mex_1p(x,F,se_opt);
%int_fcn = @(F) SE_fg_int_force_mex_1p(x,F,se_opt);

% grid
%H = cell([dim_in, 1]);
grid_t = tic;
% Grid
H = grid_fcn(f);
walltime.grid = walltime.grid + toc(grid_t);


% transform and shift, let FFT do padding
fft_t = tic;
if(se_opt.sg~=se_opt.sl)
    [G, Gres, G0]= fftnd(H, se_opt.M,se_opt.My,se_opt.Mz,se_opt.sl,...
                         se_opt.s0,se_opt.sg,se_opt.local_pad,se_opt.k0mod,3);
else
    G = fftn(H, [se_opt.M se_opt.My*se_opt.sg se_opt.Mz*se_opt.sg]);
walltime.fft = walltime.fft + toc(fft_t);
end

% scale
%scale_t = tic;
if(se_opt.sl~=se_opt.sg)
    [G, Gres, G0, scale_t] = se1p_k_scaling(G, Gres, G0, se_opt);
else
    [G, ~, ~, scale_t] = se1p_k_scaling(G,zeros(size(G)), ...
                                        zeros(se_opt.s0*se_opt.My, ...
                                              se_opt.s0*se_opt.Mz),se_opt);
end
%walltime.scale = toc(scale_t);
walltime.scale = scale_t;

% inverse shift and inverse transform
fft_t = tic;
if(se_opt.sg~=se_opt.sl)
    F = ifftnd(G,Gres,G0,se_opt.M,se_opt.My,se_opt.Mz,...
               se_opt.sl,se_opt.s0,se_opt.sg,se_opt.local_pad,se_opt.k0mod,3);
else
    F = ifftn( G ); % this should be real to eps accuracy!
    F = F(1:se_opt.Mx,1:se_opt.My,1:se_opt.M);
end
walltime.fft = walltime.fft + toc(fft_t);

u = zeros(N, 3);
int_t = tic;
u = 4*pi*int_fcn(F);
walltime.int = walltime.int + toc(int_t);

if nargout==2
    walltime.total = sum(struct2array(walltime));
    varargout{1} = walltime;
end


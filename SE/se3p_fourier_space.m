function [u varargout]= se3p_fourier_space(x, f, opt)

% initialize time array
walltime = struct('pre',0,'grid',0,'fft',0,'scale',0,'int',0);

% Check
assert(isfield(opt, 'xi'), 'xi must be given in opt struct')

% Setup vars, modify opt

se_opt = parse_params(opt);

fsize = size(f);
N = fsize(1);

% === Use vectorized code
% Gridder
pre_t = tic;
S = se3p_precomp(x,se_opt);
walltime.pre = toc(pre_t);
grid_fcn = @(f) SE_fg_grid_split_thrd_mex(x(S.perm,:),f(S.perm), ...
					 se_opt,S.zs,S.zx,S.zy,S.zz,S.idx);
% Integrator
SI = S;
iperm = @(u) u(SI.iperm,:);
int_fcn = @(F) iperm(SE_fg_int_split_mex(0,F,se_opt,... 
					SI.zs,SI.zx,SI.zy,SI.zz,SI.idx));

% === Uncomment for direct code

% grid_fcn = @(f) SE_fg_grid_mex(x,f,se_opt);
% int_fcn = @(f) SE_fg_int_mex(x,f,se_opt);

% grid
grid_t = tic;
% Grid
H = grid_fcn(f);
walltime.grid = walltime.grid + toc(grid_t);


% transform and shift, let FFT do padding
fft_t = tic;
H = fftn(H);
walltime.fft = walltime.fft + toc(fft_t);

% scale
scale_t = tic;
% k-vectors
[k1 k2 k3] = k_vec(se_opt.M,se_opt.box);
[K1 K2 K3] = ndgrid(k1,k2,k3);

% scale
k2 = K1.^2 + K2.^2 + K3.^2;
Z = exp(-(1-se_opt.eta)*k2/(4*se_opt.xi^2))./k2;
Z(1,1,1) = 0;
H = H.*Z;
walltime.scale = toc(scale_t);

% inverse shift and inverse transform
fft_t = tic;
H = ifftn( H );			% this should be real to eps accuracy!
walltime.fft = walltime.fft + toc(fft_t);

u = zeros(N, 1);
int_t = tic;
u = 4*pi*int_fcn(H);
walltime.int = walltime.int + toc(int_t);

if nargout==2
    walltime.total = sum(struct2array(walltime));
    varargout{1} = walltime;
end

% ---------------------------------------------------------------
function popt = parse_params(opt)
assert(isfield(opt,'M'));
popt = opt;
M = opt.M;
h = opt.box(1)/opt.M(1);
assert(abs( h-opt.box(2)/opt.M(2))<eps), assert(abs( h-opt.box(3)/opt.M(3))<eps)
if( isfield(opt,'P')), P = opt.P; else P = min(M); end;
if( isfield(opt,'m')), m = opt.m; else m = 1.71*sqrt(P); end;
if( isfield(opt,'w')), w = opt.w; else w = h*P/2; end;
eta = (2*w*opt.xi/m)^2;
c = 2*opt.xi^2/eta;

popt.h = h;
popt.m = m;
popt.w = w;
popt.P = P;
popt.M = M;
popt.eta = eta;
popt.c = c;

% ------------------------------------------------------------------------------
function [k1 k2 k3] = k_vec(M,box)
if (all(mod(M,2)==0))
  MM = M/2;
  k1 = (2*pi/box(1))*[0:(MM(1)-1) -MM(1):-1];
  k2 = (2*pi/box(2))*[0:(MM(2)-1) -MM(2):-1];
  k3 = (2*pi/box(3))*[0:(MM(3)-1) -MM(3):-1];

elseif(all(mod(M-1,2)==0))
  MM = (M-1)/2;
  k1 = (2*pi/box(1))*[0:MM(1) -MM(1):-1];
  k2 = (2*pi/box(2))*[0:MM(2) -MM(2):-1];
  k3 = (2*pi/box(3))*[0:MM(3) -MM(3):-1];

else error('k-vectors not computed (FIXME)');
end

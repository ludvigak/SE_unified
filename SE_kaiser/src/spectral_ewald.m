function [phi walltime phi0 walltime0] = spectral_ewald(eval_idx,x,q,opt)

fgg=true;

verb = true;
% initialize time array
walltime = struct('pre',0,'grid',0,'fft',0,'scale',0,'int',0);
walltime0= struct('pre',0,'grid',0,'fft',0,'scale',0,'int',0);

% parameters and constants
opt = parse_params(opt)

%gridding precomp
t = tic;
S = SE_FGG_precomp(x,opt.xi,opt);
grid_fcn = @(f) SE_fg_grid_split_thrd_mex(x(S.perm,:),f(S.perm),opt,...
					 S.zs,S.zx,S.zy,S.zz, S.idx);
walltime.pre = toc(t);
walltime0.pre = toc(t);

% to grid function
t=tic;
H = SE_fg_grid_kaiser_mex(x,q,opt);
walltime.grid = toc(t);

t=tic;
H0 = SE_fg_grid_mex(x,q,opt);
walltime0.grid = toc(t);

% transform and shift
t=tic;
H = fftn(H);
walltime.fft = toc(t);
t=tic;
H0 = fftn(H0);
walltime0.fft = toc(t);

% k-vectors
[k1 k2 k3] = k_vectors(opt.M,opt.box); 
[K1 K2 K3] = ndgrid(k1,k2,k3);

% scale
k2 = K1.^2 + K2.^2 + K3.^2;

t=tic;
pre = precomp(opt);
walltime.pre = walltime.pre + toc(t);
t=tic;
Z = exp(-k2/(4*opt.xi^2))./k2.*pre.F;
Z(1,1,1) = 0;
H = H.*Z;
walltime.scale = toc(t);

t=tic;
Z0 = exp(-(1-opt.eta)*k2/(4*opt.xi^2))./k2;
Z0(1,1,1) = 0;
H0 = H0.*Z0;
walltime0.scale = toc(t);


% inverse shift and inverse transform
t=tic;
H = ifftn( H );
walltime.fft = walltime.fft + toc(t);

t=tic;
H0 = ifftn( H0 );
walltime0.fft = walltime0.fft + toc(t);

% spread and integrate
t=tic;
phi = SE_fg_int_kaiser_mex(x(eval_idx,:),H,opt);
walltime.int = toc(t);

t=tic;
phi0 = SE_fg_int_mex(x(eval_idx,:),H0,opt);
walltime0.int = toc(t);

phi  = 4*pi*phi;
phi0 = 4*pi*phi0;
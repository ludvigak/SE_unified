function [phi walltime] = spectral_ewald_kaiser(eval_idx,x,q,opt)

fkg = true;

verb = true;
% initialize time array
walltime = struct('pre',0,'grid',0,'fft',0,'ifft',0,'scale',0,'int',0,'prefft',0);

% parameters and constants
opt = parse_params(opt);

if(fkg)
    pre_t = tic;
    S = SE_FKG_precomp(x,opt.xi,opt);
    grid_fcn = @(f) SE_fg_grid_split_kaiser_mex(x(S.perm,:),f(S.perm),opt,...
                                                S.zx,S.zy,S.zz, S.idx);
    walltime.pre = toc(pre_t);
    % Integrator
    SI = S;
    iperm = @(u) u(SI.iperm,:);
    int_fcn = @(F) iperm(SE_fg_int_split_kaiser_mex(0,F,opt,SI.zx, ...
                                                    SI.zy,SI.zz,SI.idx));
else
    grid_fcn = @(f) SE_fg_grid_kaiser_mex(x,f,opt);
    int_fcn = @(f) SE_fg_int_kaiser_mex(x,f,opt);
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

pre_t=tic;
pre = precomp(opt);
walltime.prefft = walltime.prefft + toc(pre_t);

Z = exp(-k2/(4*opt.xi^2))./k2.*pre.F;
Z(1,1,1) = 0;
scale_t=tic;
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

phi  = 4*pi*phi;
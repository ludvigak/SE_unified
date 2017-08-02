function [phi walltime] = se3p_fourier_space_kaiser(eval_idx,x,q,opt)

fkg = true;

verb = true;
% initialize time array
walltime = struct('pre',0,'grid',0,'fft',0,'ifft',0,'scale',0,'int',0,'prefft',0);

% parameters and constants
opt = se3p_parse_params(opt);

%gridding precomp
if(fkg)
    pre_t = tic;
    S = se3p_precomp_kaiser(x,opt.xi,opt);
    grid_fcn = @(f) SE_fg_grid_split_kaiser_mex(x(S.perm,:),f(S.perm),opt,...
                                                S.zx,S.zy,S.zz, ...
                                                S.idx);
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

% precompute kaiser window Fourier transform
pre_t=tic;
pre = se3p_window_precomp(opt);
walltime.prefft = walltime.prefft + toc(pre_t);

scale_t=tic;
H = se3p_k_scaling_kaiser(H,pre,opt);
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
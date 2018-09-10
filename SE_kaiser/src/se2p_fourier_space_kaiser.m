function [phi varargout] = se2p_fourier_space_kaiser(eval_idx,x,q,opt)

fkg = true;

verb = true;
% initialize time array
walltime = struct('pre',0,'grid',0,'fft',0,'scale',0,'int',0,'prefft',0);

% parameters and constants
se_opt = se2p_parse_params(opt);

%gridding precomp
if(fkg)
    pre_t = tic;
    S = se2p_precomp_kaiser(x,se_opt);
    grid_fcn = @(f) SE_fg_grid_split_kaiser_mex_2p(x(S.perm,:),f(S.perm),se_opt,...
                                                   S.zx,S.zy,S.zz, ...
                                                   S.idx);
    walltime.pre = toc(pre_t);
    % Integrator
    SI = S;
    iperm = @(u) u(SI.iperm,:);
    int_fcn = @(F) iperm(SE_fg_int_split_kaiser_mex_2p(0,F,se_opt,SI.zx, ...
                                                      SI.zy,SI.zz,SI.idx));
else
    grid_fcn = @(f) SE_fg_grid_kaiser_mex_2p(x,f,se_opt);
    int_fcn = @(f) SE_fg_int_kaiser_mex_2p(x,f,se_opt);
end

% to grid function
grid_t=tic;
H = grid_fcn(q);
walltime.grid = toc(grid_t);

% transform and shift
fft_t=tic;
[G Gr G0] = fftnd2p(H,se_opt);
walltime.fft = toc(fft_t);

% precompute kaiser window Fourier transform
pre_t=tic;
pre = se2p_window_precomp(se_opt);
walltime.prefft = walltime.prefft + toc(pre_t);

scale_t=tic;
[G Gr G0] = se2p_k_scaling_kaiser(G,Gr,G0,pre,se_opt);
walltime.scale = toc(scale_t);

% inverse shift and inverse transform
ifft_t=tic;
F = ifftnd2p(G,Gr,G0,se_opt);
walltime.fft = walltime.fft + toc(ifft_t);

% spread and integrate

int_t=tic;
phi = 4 * pi *int_fcn(real(F));
walltime.int = toc(int_t);

if nargout==2
    walltime.total = sum(struct2array(walltime));
    varargout{1} = walltime;
end
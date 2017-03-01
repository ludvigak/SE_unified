function [u walltime]  = SE2P_Stokes(eval_idx,x,f,xi,opt)
% SPECTRAL EWALD 2DP
% Fast Ewald method for electrostatic potential calculation, k-space part.
%
% phi = spectral_ewald_2dp(idx, x, q, xi, opt)
%
% Dag Lindbo, dag@kth.se, Jul. 2011

% initialize time array
walltime = struct('pre',0,'grid',0,'fft',0,'scale',0,'int',0);
    
% parameters and constants
opt = parse_params(opt);
[w m M Mz P grid_method spf] = unpack_params(opt);
eta = (2*w*xi/m)^2;
opt.c = 2*xi^2/eta;
opt.a = opt.wbox(3,1);
fsize = size(f);
N = fsize(1);
dim_in = fsize(2:end);

% === Use vectorized code
% Gridder
pre_t = tic;
S = SE2P_Stokes_pre(x,xi,opt);
walltime.pre = toc(pre_t);
grid_fcn = @(f) SE_fg_grid_split_thrd_mex(x(S.perm,:),f(S.perm), ...
                                          opt,S.zs,S.zx,S.zy,S.zz,S.idx);

% Integrator
SI = S;
iperm = @(u) u(SI.iperm,:);
int_fcn = @(F) iperm(SE_fg_int_split_mex(x(SI.iperm,:),F, ...
                                         opt,SI.zs,SI.zx,SI.zy,SI.zz,SI.idx));

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
    H{i} = fftshift( fftn(tmp, [M M Mz*spf]) );
    walltime.fft = walltime.fft + toc(fft_t);
end

% scale
scale_t = tic;
[G1 G2 G3] = mul_k_space_kernel_mex(H{1}, H{2}, H{3}, opt, xi, eta);
walltime.scale = toc(scale_t);
G{1} = G1;
G{2} = G2;
G{3} = G3;
dim_out = numel(G);

% inverse shift and inverse transform
for i=1:dim_out
    fft_t = tic;
    G{i} = ifftn( ifftshift(G{i}));
    walltime.fft = walltime.fft + toc(fft_t);
end

% integrate
u = zeros(N, dim_out);
for i=1:dim_out
    int_t = tic;
    u(:,i) = int_fcn(G{i}(:, :, 1:Mz));
    walltime.int = walltime.int + toc(int_t);
end
if nargout==2
    walltime.total = sum(struct2array(walltime));
    varargout{1} = walltime;
end

% ------------------------------------------------------------------------------
function [G1 G2 G3] = mul_k_space_kernel(g1, g2, g3, opt, xi, eta, op_B)
k1    = k_vectors(opt.M, opt.L, 1);
k2    = k_vectors(opt.M, opt.L, 1);
kappa = k_vectors(opt.Mz,opt.Lz,opt.sampling_factor);
[K1 K2 KAPPA] = ndgrid(k1,k2,kappa);

G1 = zeros(size(g1));
G2 = zeros(size(g2));
G3 = zeros(size(g3));
for i1 = 1:length(k1)
    for i2 = 1:length(k2)
        for i3 = 1:length(kappa)
            k = [K1(i1,i2,i3) K2(i1,i2,i3) KAPPA(i1,i2,i3)];
            if all(k(1:2) == 0)
                G1(i1,i2,i3) = 0; G2(i1,i2,i3) = 0; G3(i1,i2,i3) = 0;
                continue,
            end
            G = exp(-(1-eta)*norm(k,2)^2/(4*xi^2))*op_B(k,xi)*...
                [g1(i1,i2,i3); g2(i1,i2,i3); g3(i1,i2,i3)];
            G1(i1,i2,i3) = G(1); G2(i1,i2,i3) = G(2); G3(i1,i2,i3) = G(3);
        end
    end
end

% ------------------------------------------------------------------------------
function [w m M Mz P grid_method spf] = unpack_params(opt)
w = opt.w;
m = opt.m;
M = opt.M;
Mz = opt.Mz;
P = opt.P;
grid_method = opt.grid_method;
spf = opt.sampling_factor;

% ------------------------------------------------------------------------------
function [ks idxz] = k_vectors(M,L,q)

% q: oversampling factor (1 = no oversampling)
Mq = q*M;

if mod(Mq,2)==0
    k = (-Mq/2):(Mq/2-1);
    idxz = Mq/2+1;
else
    k = -(Mq-1)/2:(Mq-1)/2;
    idxz = (Mq-1)/2+1;
end
ks = 2*pi*k/(L*q);
assert(abs(k(idxz))<eps)

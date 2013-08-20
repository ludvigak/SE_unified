function [u stats]  = SE2P_Stokes(eval_idx,x,f,xi,opt,varargin)
% SPECTRAL EWALD 2DP
% Fast Ewald method for electrostatic potential calculation, k-space part.
%
% phi = spectral_ewald_2dp(idx, x, q, xi, opt)
%
% Dag Lindbo, dag@kth.se, Jul. 2011

verb = false;

% parameters and constants
opt = parse_params(opt);
[w m M Mz P grid_method spf] = unpack_params(opt);
eta = (2*w*xi/m)^2;
opt.c = 2*xi^2/eta;
opt.a = opt.wbox(3,1);

if nargin == 5
    static_fgg=false;
elseif nargin == 6
    static_fgg=true;
    sdat=varargin{1};
end

cprintf(verb,'[ 2DP SPEC EWALD ] ')
cprintf(verb, 'Grid: [%d %d %d(%d)]\tGaussian: {P = %d, w=%.2f, m=%.1f}\n',...
        M, M, Mz, spf, P, w, m)
cprintf(verb,'[ 2DP SPEC EWALD ] Gridding method: %s {eta=%.2f, c=%.2f}\n', ...
        grid_method, eta, opt.c)

% to grid function
tic
if(static_fgg)

    x = x(sdat.perm,:);
    f = f(sdat.perm,:);
    
    H1 = SE_fg_grid_split_mex(x,f(:,1),opt,...
                              sdat.zs,sdat.zx,sdat.zy,sdat.zz,sdat.idx);
    H2 = SE_fg_grid_split_mex(x,f(:,2),opt,...
                              sdat.zs,sdat.zx,sdat.zy,sdat.zz,sdat.idx);
    H3 = SE_fg_grid_split_mex(x,f(:,3),opt,...
                              sdat.zs,sdat.zx,sdat.zy,sdat.zz,sdat.idx);
else
    H1 = SE_fg_grid_mex(x,f(:,1),opt);
    H2 = SE_fg_grid_mex(x,f(:,2),opt);
    H3 = SE_fg_grid_mex(x,f(:,3),opt);
end
stats.wtime_grid = toc();

% transform and shift
tic
G1 = fftshift( fftn(H1, [M M Mz*spf]) );
G2 = fftshift( fftn(H2, [M M Mz*spf]) );
G3 = fftshift( fftn(H3, [M M Mz*spf]) );
stats.wtime_fft = toc();

% multiply with k-space kernel
tic
[G1 G2 G3] = mul_k_space_kernel_mex(G1, G2, G3, opt, xi, eta);
stats.wtime_poisson = toc();

% inverse shift and inverse transform
tic
F1 = ifftn( ifftshift( G1 )); 
F2 = ifftn( ifftshift( G2 )); 
F3 = ifftn( ifftshift( G3 )); 

F1 = F1(:,:,1:Mz);
F2 = F2(:,:,1:Mz);
F3 = F3(:,:,1:Mz);
stats.wtime_fft = stats.wtime_fft + toc();

% integrate
tic
u = zeros(length(eval_idx),3);
if(static_fgg)
    u(:,1) = SE_fg_int_split_mex(x,F1,opt,...
                                 sdat.zs,sdat.zx,sdat.zy,sdat.zz,sdat.idx);
    u(:,2) = SE_fg_int_split_mex(x,F2,opt,...
                                 sdat.zs,sdat.zx,sdat.zy,sdat.zz,sdat.idx);
    u(:,3) = SE_fg_int_split_mex(x,F3,opt,...
                                 sdat.zs,sdat.zx,sdat.zy,sdat.zz,sdat.idx);
    u = u(sdat.iperm,:);
else
    u(:,1) = SE_fg_int_mex(x(eval_idx,:),F1,opt);
    u(:,2) = SE_fg_int_mex(x(eval_idx,:),F2,opt);
    u(:,3) = SE_fg_int_mex(x(eval_idx,:),F3,opt);
end
stats.wtime_int = toc();
stats.wtime_total = stats.wtime_grid + stats.wtime_fft + ...
    stats.wtime_poisson + stats.wtime_int;
stats.wtime=[stats.wtime_grid stats.wtime_fft ...
    stats.wtime_poisson stats.wtime_int];

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

% ------------------------------------------------------------------------------
function st = make_stats(opt,wtime)
st.wall_time_breakdown = wtime;
st.wall_time_labels = {'grid','FFT','Poisson','IFFT','int'};
st.grid_size = [opt.M opt.M opt.Mz];
st.grid_prod = prod(st.grid_size);
st.fft_size = [opt.M opt.M opt.Mz*opt.sampling_factor];
st.fft_prod = prod(st.fft_size);


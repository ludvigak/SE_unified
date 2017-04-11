function phi  = spectral_ewald_2P(eval_idx,x,q,xi,opt,varargin)
% SPECTRAL EWALD 2DP
% Fast Ewald method for electrostatic potential calculation, k-space part.
%
% phi = spectral_ewald_2dp(idx, x, q, xi, opt)
%
% Dag Lindbo, dag@kth.se, Dec. 2010

verb = true;

% parameters and constants
opt = parse_params(opt);
[w m M Mz P grid_method spf] = unpack_params(opt);
eta = (2*w*xi/m)^2;
opt.c = 2*xi^2/eta;
opt.a = opt.wbox(3,1);

cprintf(verb,'[ 2DP SPEC EWALD ] ')
cprintf(verb, 'Grid: [%d %d %d(%d)]\tGaussian: {P = %d, w=%.2f, m=%.1f}\n',...
        M, M, Mz, spf, P, w, m)
cprintf(verb,'[ 2DP SPEC EWALD ] Gridding method: %s {eta=%.2f, c=%.2f}\n', ...
        grid_method, eta, opt.c)


if nargin == 5
    static_fgg=false;
elseif nargin == 6
    static_fgg=true;
    sdat=varargin{1};
end

% to grid function
if(static_fgg)
    x = x(sdat.perm,:);
    q = q(sdat.perm);
    H = SE_fg_grid_split_mex_2p(x,q,opt,...
        sdat.zs,sdat.zx,sdat.zy,sdat.zz,sdat.idx);
else
    H  = SE_fg_grid_mex_2p(x,q,opt);
end


% transform and shift
G = fftn(H, [M M Mz*spf]);

% k-vectors
[k1 zidx1] = k_vectors(M, opt.L, 1);
[k2 zidx2] = k_vectors(M, opt.L, 1);
kappa      = k_vectors(Mz,opt.Lz,spf);
[K1 K2 KAPPA] = ndgrid(k1,k2,kappa);

% scale
k2 = K1.^2 + K2.^2 + KAPPA.^2;
Z = exp(-(1-eta)*k2/(4*xi^2))./k2;
Z(zidx1,zidx2,:) = 0;
G = G.*Z;

% inverse shift and inverse transform
F = ifftn( G ); % this should be real to eps accuracy!
F = F(:,:,1:Mz);

% integrate
if(static_fgg)
    assert(length(eval_idx)==size(x,1),'fixme')
    phi = SE_fg_int_split_mex_2p(x,F,opt,...
                                 sdat.zs,sdat.zx,sdat.zy,sdat.zz,sdat.idx);
    phi = phi(sdat.iperm,:);
else
    phi = SE_fg_int_mex_2p(x(eval_idx,:),F,opt);
    
end
phi = 4*pi*phi/(opt.L^2);

% done

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
else
  k = -(Mq-1)/2:(Mq-1)/2;
end
k = fftshift(k); % standard reordering
idxz = 1;
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


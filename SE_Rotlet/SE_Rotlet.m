function varargout = SE_Rotlet(xe,x,f,xi,opt)
% u = SE_Rotlet(x_targets, x_sources, f_sources, xi, opt)
%   Returns the rotlet potential
%
% [U1, U2, U3] = SE_Rotlet(x_targets, x_sources, f_sources, xi, opt)
%   Returns the grid with Fourier coefficients, 
%   :math:`U_i\in\mathbb{C}^{M_1\times M_2 \times M_3}`
%
% To compute u from U1,U2,U3:
%
% .. code:: matlab
%
%  F{1} = real( ifftn( ifftshift( U1 )));
%  F{2} = real( ifftn( ifftshift( U2 )));
%  F{3} = real( ifftn( ifftshift( U3 )));
%  u = SE_fg_int(xe, F, opt);
% 
% opt must be identical to that passed to SE_Rotlet
%
% TODO: Add possibility of passing precomputed structures

verb = false;

% parameters and constants
opt.xi = xi;
opt = parse_params(opt);
N = size(x, 1);

% Use vectorized code
% Gridder
S = SE_FGG_precomp(x,xi,opt);
grid_fcn = @(f) SE_fg_grid_split_mex(x(S.perm,:),f(S.perm),opt,S.zs,S.zx,S.zy,S.zz, ...
                                     S.idx);
% Integrator
SI = SE_FGG_precomp(xe,xi,opt);
iperm = @(u) u(SI.iperm,:);
int_fcn = @(F) iperm(SE_fg_int_split_mex(0,F,opt,SI.zs,SI.zx,SI.zy,SI.zz,SI.idx));

% Uncomment for direct code
%grid_fcn = @(f) SE_fg_grid_mex(x,f,opt);
%int_fcn = @(f) SE_fg_int_mex(x(eval_idx,:),f,opt);

% to grid function
H1 = grid_fcn(f(:,1));
H2 = grid_fcn(f(:,2));
H3 = grid_fcn(f(:,3));

% transform and shift
H1 = fftshift( fftn(H1) );
H2 = fftshift( fftn(H2) );
H3 = fftshift( fftn(H3) );

% Scale
[k1 k2 k3] = k_vectors(opt.M, opt.box); 
[K1 K2 K3] = ndgrid(k1,k2,k3);
Ksq = K1.^2 + K2.^2 + K3.^2;
B = exp(-(1-opt.eta)*Ksq/(4*xi^2))./Ksq; % keep it real
B(k1==0, k2==0, k3==0) = 0;
% G = B*(HxK)
G1 = -1i*4*pi*B.*(H2.*K3 - H3.*K2);
G2 = -1i*4*pi*B.*(H3.*K1 - H1.*K3);
G3 = -1i*4*pi*B.*(H1.*K2 - H2.*K1);

if nargout > 1
    varargout{1} = G1;
    varargout{2} = G2;
    varargout{3} = G3;
    return
end

% inverse shift and inverse transform
F1 = real( ifftn( ifftshift( G1 )));
F2 = real( ifftn( ifftshift( G2 )));
F3 = real( ifftn( ifftshift( G3 )));

u = zeros(size(xe));
u(:,1) = int_fcn(F1);
u(:,2) = int_fcn(F2);
u(:,3) = int_fcn(F3);

varargout{1} = u;

% ------------------------------------------------------------------------------
function [k] = k_vec(M,L)
if mod(M,2)==0
    MM = M/2;
    k = (2*pi/L)*[-MM:(MM-1)];
elseif mod(M-1,2)==0
    MM = (M-1)/2;
    k = (2*pi/L)*[-MM:MM];
else error('k-vectors not computed');
end

% ------------------------------------------------------------------------------
function [k1 k2 k3] = k_vectors(M,box)
k1 = k_vec(M(1), box(1));
k2 = k_vec(M(2), box(2));
k3 = k_vec(M(3), box(3));


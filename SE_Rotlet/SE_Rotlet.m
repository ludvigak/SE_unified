function u = SE_Rotlet(xe,x,f,xi,opt)
% u = SE_Rotlet(x_targets, x_sources, f_sources, xi, opt)

% TODO: Add possibility of passing precomputed structures

verb = false;

% parameters and constants
opt = parse_params(opt);
[w m M P] = unpack_params(opt);
eta = (2*w*xi/m)^2;
opt.c = 2*xi^2/eta;
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

cprintf(verb, 'M = [%d %d %d] P = %d m=%d w=%f\n',M,P,m,w);
cprintf(verb, 'eta = %f\t a=%f\n', eta, pi^2/opt.c);

% transform and shift
H1 = fftshift( fftn(H1) );
H2 = fftshift( fftn(H2) );
H3 = fftshift( fftn(H3) );

% Scale
[k1 k2 k3] = k_vectors(M,opt.box); 
[K1 K2 K3] = ndgrid(k1,k2,k3);
Ksq = K1.^2 + K2.^2 + K3.^2;
B = exp(-(1-eta)*Ksq/(4*xi^2))./Ksq; % keep it real
B(k1==0, k2==0, k3==0) = 0;
% G = B*(HxK)
G1 = -1i*4*pi*B.*(H2.*K3 - H3.*K2);
G2 = -1i*4*pi*B.*(H3.*K1 - H1.*K3);
G3 = -1i*4*pi*B.*(H1.*K2 - H2.*K1);


% inverse shift and inverse transform
F1 = real( ifftn( ifftshift( G1 )));
F2 = real( ifftn( ifftshift( G2 )));
F3 = real( ifftn( ifftshift( G3 )));

u = zeros(size(xe));
u(:,1) = int_fcn(F1);
u(:,2) = int_fcn(F2);
u(:,3) = int_fcn(F3);


% ------------------------------------------------------------------------------
function p = parse_params(opt)

% check that we have all mandatory options
assert(isfield(opt,'M'))
assert(isfield(opt,'P'))
assert(isfield(opt,'box'))

% verify all assumptions on parameters

% step size
L = opt.box(1);
h = L/opt.M(1);
assert(abs(opt.box(2)/opt.M(2) - opt.box(1)/opt.M(1)) < eps)
assert(abs(opt.box(3)/opt.M(3) - opt.box(2)/opt.M(2)) < eps)

% Gaussian
P = opt.P;
if( isfield(opt,'m')), m = opt.m; else m = 0.9*sqrt(pi*P); end;
w = h*(P-1)/2;

% collect
p.M=opt.M;
p.P = P;
p.w = w;
p.m = m;
p.box = opt.box;
p.L = L;
p.h = h;


% ------------------------------------------------------------------------------
function [w m M P] = unpack_params(opt)
w = opt.w;
m = opt.m;
M = opt.M;
P = opt.P;

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


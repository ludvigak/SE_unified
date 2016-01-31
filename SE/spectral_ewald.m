function phi  = spectral_ewald(eval_idx,x,q,xi,opt)
% SPECTRAL EWALD
% Fast Ewald method for electrostatic potential calculation, k-space part.
%
% phi = spectral_ewald(idx, x, q, xi, box, opt)
%
% :return: **phi** -- k-space part of periodic potential
%
% :param idx: evaluation indices (particle indices)
% :param x:   particle positions (N-by-3)
% :param q:   particle charge (N-by-1)
% :param xi:  Ewald parameter
% :param opt: Options:
% :param opt.M:           grid size, eg. [31 31 31] (required)
% :param opt.P:           Gaussian support    (required)
% :param opt.box:         Periodic cell, eg. [1 1 1]
% :param opt.m:           Gaussian shape function (default: m=m(P))
% :param opt.w:           Gaussian width (default: w=h*P/2)
%

verb = true;

% parameters and constants
[M P m w grid_method] = parse_params(opt);
eta = (2*w*xi/m)^2;
c = 2*xi^2/eta;
fgg_opt.c = c; fgg_opt.M=M; fgg_opt.P=P; fgg_opt.box=opt.box;

cprintf(verb, '[ SE ] Grid: [%d %d %d]\tGaussian: {P = %d, w=%.2f, m=%d}\n',...
    M, P, w, m)
cprintf(verb, '[ SE ] Gridding method: %s {eta=%.2f, c=%.2f}\n', ...
    grid_method, eta, c)

% sorting by closest grid index
assert(length(eval_idx)==length(q))
[x perm] = sort_by_grid_index(x,opt.box(1)/M(1),M(1));
iperm(perm)=1:length(perm);
q = q(perm);

% to grid function
H = SE_fg_grid_mex(x,q,fgg_opt);

% transform and shift
H = fftn(H);

% k-vectors
[k1 k2 k3] = k_vec(M,opt.box); 
[K1 K2 K3] = ndgrid(k1,k2,k3);

% scale
k2 = K1.^2 + K2.^2 + K3.^2;
Z = exp(-(1-eta)*k2/(4*xi^2))./k2;
Z(1,1,1) = 0;
H = H.*Z;

% inverse shift and inverse transform
H = ifftn( H );

% spread and integrate
phi = SE_fg_int_mex(x(eval_idx,:),H,fgg_opt);

phi = 4*pi*phi;
phi = phi(iperm); % unsorting

% ------------------------------------------------------------------------------
function [M P m w meth] = parse_params(opt)
assert(isfield(opt,'M'));

M = opt.M;
h = opt.box(1)/opt.M(1);
assert(abs( h-opt.box(2)/opt.M(2))<eps), assert(abs( h-opt.box(3)/opt.M(3))<eps)
if( isfield(opt,'P')), P = opt.P; else P = min(M); end;
if( isfield(opt,'m')), m = opt.m; else m = 1.71*sqrt(P); end;
if( isfield(opt,'w')), w = opt.w; else w = h*P/2; end;
if( isfield(opt,'grid_method')), meth = opt.grid_method; 
else meth = 'fgg mex'; end;

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

function [xs p] = sort_by_grid_index(x,h,M)

sub = zeros(size(x));
sub(:,1) =  mod(round(x(:,1)/h),M)+1;
sub(:,2) =  mod(round(x(:,2)/h),M)+1;
sub(:,3) =  mod(round(x(:,3)/h),M)+1;

idx = sub(:,1)*M*M + sub(:,2)*M + sub(:,3);

[~, p] = sort(idx);
xs = x(p,:);

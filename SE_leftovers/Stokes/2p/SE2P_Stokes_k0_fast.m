% ------------------------------------------------------------------------------
% Ewadl sum under 2dp conditions, singular k-space part
function u = SE2P_Stokes_k0_fast(idx, x, f, opt)

% parameters
M = opt.k0_M;
xi = opt.xi;
box = opt.box;
if isfield(opt,'verb'), verbose = opt.verb; else verbose = true; end
if isfield(opt,'method'), meth = ['_' opt.method]; else meth = '_mex'; end
if isfield(opt,'interp_meth'), interp_method = opt.interp_meth; 
else interp_method = 'spectral'; end

if verbose
    fprintf('[SE2P STOKES (K0F)] N=(%d,%d) M=%d, Method: %s, %s\n',...
        size(x,1),length(idx),M,interp_method,meth);
end

[zm I]= sort(x(idx,3)); % observation points sorted

switch interp_method
 case 'spectral'
    z = gauss_points(0,box(3),M); % gauss/lobotto points
    %H1 = SE2P_k0_kernel_mex(x,f(:,1),xi,z);
    %H2 = SE2P_k0_kernel_mex(x,f(:,2),xi,z); 
    [H1 H2] = SE2P_k0_kernel_mex(x,f(:,1),f(:,2),xi,z);
    u1 = chebinterp1(H1,zm,0,box(3));
    u2 = chebinterp1(H2,zm,0,box(3));
  otherwise
    error('invalid interp method')
end

IV(I) = 1:length(I); % inverse permutation to undo the sorting
c = -4*sqrt(pi)/(opt.box(1)*opt.box(2));
u = c*[u1(IV) u2(IV) zeros(size(u1))];

%phi = -2*phi(IV)*sqrt(pi)/(box(1)*box(2));
%done

function H = SE2P_k0_kernel(x,q,xi,z)
% Evaluate sum over particles with respect to (grid-) points z.
% Has MEX equivalent;
H = zeros(size(z));

for j = 1:length(z)                 % grid points
    phi_j = 0;
    for n = 1:size(x,1)               % particles
        zjn = x(n,3)-z(j);
        phi_j = phi_j + q(n)*(exp(-xi^2*zjn^2)/(2*xi) + ...
                              sqrt(pi)*zjn*erf(xi*zjn));
    end
    H(j) = phi_j;
end

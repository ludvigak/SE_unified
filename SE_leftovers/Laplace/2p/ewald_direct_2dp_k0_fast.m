% ------------------------------------------------------------------------------
% Ewadl sum under 2dp conditions, singular k-space part
function phi = ewald_direct_2dp_k0_fast(idx, x, q, opt)

% parameters
M = opt.k0_M;
xi = opt.xi;
box = opt.box;
if isfield(opt,'verb'), verbose = opt.verb; else verbose = true; end
if isfield(opt,'method'), meth = ['_' opt.method]; else meth = ''; end
if isfield(opt,'interp_meth'), interp_method = opt.interp_meth; 
else interp_method = 'spectral'; end

if verbose
    fprintf('[2DP EWALD (K0F )] N=(%d,%d) M=%d, Method: %s, %s\n',...
        size(x,1),length(idx),M,interp_method,meth);
end

[zm I]= sort(x(idx,3)); % observation points sorted

switch interp_method
 case 'spectral'
    z = gauss_points(0,box(3),M); % gauss/lobotto points
    eval(['H = ewald_direct_2dp_k0_kernel' meth '(x,q,xi,z);']);
    phi = chebinterp1(H,zm,0,box(3));
  
 case 'spline'
    z = linspace(0,box(3),M)';  % equidistant grid
    eval(['H = ewald_direct_2dp_k0_kernel' meth '(x,q,xi,z);']);
    phi = interp1(z,H,zm,'spline');

end

IV(I) = 1:length(I); % inverse permutation to undo the sorting
phi = -2*phi(IV)*sqrt(pi)/(box(1)*box(2));
%done

function H = ewald_direct_2dp_k0_kernel(x,q,xi,z)
% Evaluate sum over particles with respect to (grid-) points z.
% Has MEX equivalent;
H = zeros(size(z));

for j = 1:length(z)                 % grid points
    phi_j = 0;
    for n = 1:size(x,1)               % particles
        zjn = x(n,3)-z(j);
        phi_j = phi_j + q(n)*(exp(-xi^2*zjn^2)/xi + sqrt(pi)*zjn*erf(xi*zjn));
    end
    H(j) = phi_j;
end

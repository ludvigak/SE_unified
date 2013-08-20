% ------------------------------------------------------------------------------
% Ewadl sum under 2dp conditions, singular k-space part
function u = SE2P_Stokes_k0_fast_apply(idx, x, f, G, opt)

% parameters
M = opt.k0_M;
box = opt.box;
if isfield(opt,'verb'), verbose = opt.verb; else verbose = true; end

cprintf(verbose,'[SE2P STOKES (K0F)] Apply N=(%d,%d) M=%d... ',...
        size(x,1),length(idx),M);

[zm I]= sort(x(idx,3)); % observation points sorted

% sum over f
H1 = G*f(:,1);
H2 = G*f(:,2);

% interpolate
u1 = chebinterp1(H1,zm,0,box(3));
u2 = chebinterp1(H2,zm,0,box(3));

IV(I) = 1:length(I); % inverse permutation to undo the sorting
c = -4*sqrt(pi)/(opt.box(1)*opt.box(2));
u = c*[u1(IV) u2(IV) zeros(size(u1))];

cprintf(verbose, 'Done \n')
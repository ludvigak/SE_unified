function phi = stresslet_direct_fd_zero( idx, xvec, fvec, nvec, L)
% Component for :math:`k=0` of stresslet Fourier sum.
%    
% :param idx: index of target positions (Nx1)
% :param xvec: source positions (Nx3)
% :param fvec: source strengths (Nx3)
% :param nvec: normal vectors (Nx3)
% :param L: box size (1x3)
% :returns: **phi_k0**  -- (Nx3)


VERBOSE = 0;

cprintf(VERBOSE,'\tComputing frequency domain sum (DIRECT SUMMATION). Zero mode.\n');

% We have r = x - xs,
% compute sums for x and xs separately

% Create lists for point-wise dot products
% f = f(xs)
% n = n(xs)
fdotxs = sum(fvec.*xvec,2); % f_j xs_j
ndotxs = sum(nvec.*xvec,2); % n_j xs_j
ndotf = sum(nvec.*fvec,2);  % n_j f_j

V = prod(L);

tmp = sum( bsxfun(@times,xvec,ndotf),1); % xs_i n_j  f_j
phi = (8*pi/V)*repmat( tmp, numel(idx), 1 );

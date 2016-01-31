function phi = stresslet_direct_fd_zero( idx, xvec, fvec, nvec, L)
% k=0 part of Stresslet k-space summation.
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

% TEST
% Check if we can get it to work using a different k=0, based on mean flow
tmp = sum( bsxfun(@times,xvec,ndotf),1); % xs_i n_j  f_j
phi = (8*pi/V)*repmat( tmp, numel(idx), 1 );

return
% END TEST

% Constant term
tmpxs = sum( ...
    bsxfun(@times,nvec,fdotxs) + ... n_i f_j xs_j
    bsxfun(@times,fvec,ndotxs) + ... f_i n_j xs_j
    bsxfun(@times,xvec,ndotf) ...   xs_i n_j  f_j
    ,1);

% Sums for point-dependent term, create matrix A
% nf_ji = sum n_j f_i
nf = zeros(3,3);
for j=1:3
    for m=1:3
        nf(j,m) = sum(nvec(:,j).*fvec(:,m));
    end
end
A = nf + nf' + sum(ndotf)*eye(3);
% x*A gives:
% x_j nf_ji  = x_j sum n_j f_i
% x_j nf'_ji = x_j sum f_j n_i
% x_j d_ji sum n_k f_k = x_i f_j n_j

% For stresslet distribution on closed surface,
% it turns out that A=0 (numerically)
% UPDATE: THIS IS WRONG!

% OK, so we say we are only interested in the periodic part,
% so skipping A
if true
%     warning('A=0');
    A=A*0;
end

phi = xvec(idx,:)*A;
phi = bsxfun(@minus, phi, tmpxs);
phi = -8*pi/(5*V)*phi;

if numel(idx)>1 && ~(all(all(diff(phi)==0)))
    error('k=0 component not constant')
end
% --------------------------------------------------
% Non-vectorized code (readable, but N^2 slow)
cmp_old = 0;
if cmp_old
    phi_ref = phi;
    noeval=length(idx);
    phi=zeros(noeval,3);
    % for targets
    for ii=1:noeval
        is=idx(ii);
        tmp = [0 0 0];
        % for sources
        for it = 1:size(xvec,1)
            r=xvec(is,:)-xvec(it,:);
            f=fvec(it,:);
            n=nvec(it,:);

            fr = f(1)*r(1) + f(2)*r(2) + f(3)*r(3);
            nr = n(1)*r(1) + n(2)*r(2) + n(3)*r(3);
            nf = n(1)*f(1) + n(2)*f(2) + n(3)*f(3);

            tmp = tmp + n*fr + f*nr + r*nf;        
        end
        phi(ii,:) = phi(ii,:) + tmp;
    end;
    V = prod(L);
    phi = -8*pi/(5*V)*phi;

    norm(phi(:)-phi_ref(:), inf) / norm(phi(:),inf)
end

    
end


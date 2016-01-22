function phi = rotlet_direct_fd(xe, x, f, xi, L, kmax)
% Ewald summation for Stokes rotlet --  k-space part.
%
% phi = rotlet_direct_fd(xe, x, f, xi, L, kmax)
%        Evaluate potential, phi, at points x(idx,:). 
%        xe   --  target positions    M-by-3
%        x    --  source positions    N-by-3
%        f    --  source strengths    N-by-3
%        xi   --  Ewald parameter
%        L    --  box size            1-by-3
%        kmax --  k-space truncation
%

k = 2*pi*(-kmax:kmax); %modes
[k1 k2 k3] = ndgrid(k/L(1),k/L(2),k/L(3));
k = [k1(:) k2(:) k3(:)];
k(sum(abs(k),2) == 0,:) = [];
Nk = size(k,1);

%fprintf('\tComputing frequency domain sum (DIRECT SUMMATION). Modes: %d\n', Nk);
%disp(['WITH xi=' num2str(xi) '.']);

noeval=size(xe,1);;  
phi=zeros(noeval,3); 
parfor j = 1:Nk
    kj = k(j,:);
    k2 = kj(1)^2+kj(2)^2+kj(3)^2;
    gaussian = exp(-k2/(4*xi^2));            
    Bi = 4*pi*kj/k2;
    tmp = zeros(noeval,3);
    for ii=1:noeval
        r = bsxfun(@minus, xe(ii,:), x);
        kr = sum(bsxfun(@times, kj, r), 2);
        nsum = sum(bsxfun(@times, f, exp(-1i*kr)), 1);
        crossprod = real(1i * cross(nsum, Bi));    
        tmp(ii,:) = crossprod*gaussian;
    end
    phi = phi + tmp;
end;
phi = phi/prod(L);

end
function phi = rotlet_direct_fd( idx, x, f, xi, L, kmax)
% Ewald summation for Stokes rotlet --  k-space part.
%
% phi = rotlet_direct_fd(idx, x, f, xi, L, kmax)
%        Evaluate potential, phi, at points x(idx,:). 
%        x    --  positions           N-by-3
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

noeval=length(idx);  
phi=zeros(noeval,3); 
for ii=1:noeval
  m=idx(ii); 
  tmp = [0 0 0];
  parfor j = 1:Nk
    z = 0;
    kj = k(j,:);
    k2 = kj(1)^2+kj(2)^2+kj(3)^2;
    Bi = 4*pi*kj/k2;
    r = bsxfun(@minus, x(m,:), x);
    kr = sum(bsxfun(@times, kj, r), 2);
    nsum = sum(bsxfun(@times, f, exp(-1i*kr)), 1);
    crossprod = real(1i * cross(nsum, Bi));    
    gaussian = exp(-k2/(4*xi^2));
    tmp = tmp + crossprod*gaussian;
    %%    cprintf(mod(j,1000)==0,'\t%3d \t %5.2f%% \r',j, 100*j/Nk);
  end
  phi(ii,:) = phi(ii,:) + tmp;
end;
phi = phi/prod(L);
%disp(['phi = ']); 
%disp(phi)

fprintf('\n')
end
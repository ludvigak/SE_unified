function phi = stresslet_direct_fd( idx, x, f, nvec, xi, L, kmax)
% Ewald summation for Stokes --  k-space part.
%
% phi = stresslet_direct_fd( m, x, f, nvec, xi, L, kmax)
%        Evaluate potential, phi, at points x(idx). 
%        x    --  positions           N-by-3
%        f    --  source strengths    N-by-3
%        nvec --  normal vectors   N-by-3
%        xi   --  Ewald parameter
%        L    --  box size            1-by-3
%        kmax --  k-space truncation
%
% Example:
%   [x f] = generate_state(100,[1 1 1]);
%   stresslet_direct_fd(2,x,f,nvec,xi, [1 1 1], 10)
%

k = 2*pi*(-kmax:kmax); %modes
[k1 k2 k3] = ndgrid(k/L(1),k/L(2),k/L(3));
k = [k1(:) k2(:) k3(:)];
k(sum(abs(k),2) == 0,:) = [];
Nk = size(k,1);

fprintf('\tComputing frequency domain sum (DIRECT SUMMATION). Modes: %d\n', Nk);
disp(['WITH xi=' num2str(xi) '.']);

noeval=length(idx);  
phi=zeros(noeval,3); 
for ii=1:noeval
  m=idx(ii); 
  tmp = [0 0 0];
  parfor j = 1:Nk
    z = 0;
    Bi = stresslet_op_fd( k(j,:), xi);
    for n = 1:size(x,1)
        
      % Old notation
      %z = z + sin(dot(k(j,:), x(n,:) - x(m,:)))*stresslet_op_fd_old(k(j,:), nvec(n,:),xi)*f(n,:)';
      
      r=x(m,:)-x(n,:);
      % Full notation
      z = z + real( 1i*...
            [	nvec(n,:)*Bi(:,:,1)*f(n,:)'
                nvec(n,:)*Bi(:,:,2)*f(n,:)'
                nvec(n,:)*Bi(:,:,3)*f(n,:)' ]*...
            exp(-1i*dot(k(j,:), r)) );
      
% % %       Real part notation
%       z = z + ...
%             [	nvec(n,:)*Bi(:,:,1)*f(n,:)'
%                 nvec(n,:)*Bi(:,:,2)*f(n,:)'
%                 nvec(n,:)*Bi(:,:,3)*f(n,:)' ]*...
%             sin(dot(k(j,:), r));

    end
    gaussian = exp(-sum(k(j,:).^2)/(4*xi^2));
    tmp = tmp + gaussian*z';

    %%    cprintf(mod(j,1000)==0,'\t%3d \t %5.2f%% \r',j, 100*j/Nk);
  end
  phi(ii,:) = phi(ii,:) + tmp;
end;
phi = phi/prod(L);
%disp(['phi = ']); 
%disp(phi)

fprintf('\n')
end

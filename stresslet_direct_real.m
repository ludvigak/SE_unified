function [phi shellnorms] = stresslet_direct_real( idx, x, f, nvec, xi, L, nbox, TOL)
% Ewald summation for the stresslet -- Real space part.
%
% phi = stokes_ewald_direct_real( m, x, f, xi, op_A, nbox)
%        Evaluate potential, phi, at points x(idx). 
%        x    --  positions          N-by-3
%        nvec --  normal vectors   N-by-3
%        f    --  source strengths   N-by-3
%        xi   --  Ewald parameter
%        nbox --  periodic repications
%
% Example:
%   [x f] = generate_state(100,[1 1 1]);
%   stresslet_direct_real(1, x, f, nvec, 2, [1 1 1], 10)

p = -nbox:nbox; %boxes
[p1 p2 p3] = ndgrid(p,p,p);
p = [p1(:) p2(:) p3(:)];
Np = size(p,1);

pshell = max(abs(p'));
[pshell, I] = sort(pshell);
pshell = pshell';
p = p(I,:)*diag(L);

shellnorms=[];

fprintf('\tComputing real space sum. Periodic images: %d\n', Np);
noeval=length(idx);  
phi=zeros(noeval,3); 
for ii=1:noeval
  m=idx(ii); 
  converged=0;
  for shell_no=0:nbox    
    tmp = [0 0 0];
    shell_indices = find(pshell==shell_no);
    parfor j = shell_indices'
%     for j = 1:Np
      for n = 1:size(x,1) % particles
        if all(p(j,:)==0) && n==m % remove self interaction
          continue
        end
        r=x(m,:)-x(n,:);
        tmp = tmp + (stresslet_op_real( r + p(j,:), nvec(n,:), xi)*f(n,:)')';
      end
    %    cprintf(mod(j,1000)==0,'\t%3d \t %5.2f%% \r',j, 100*j/Np);
    end
    phi(ii,:) = phi(ii,:) + tmp;
    shellnorm = norm(tmp);
    shellnorms(end+1)=shellnorm;
    fprintf('Shell %d: contrib=%g\n',shell_no,shellnorm);
    if shellnorm<TOL && shell_no>3
        converged=1;
        break
    end
  end
  if converged==1
      fprintf('Direct sum converged TOL=%g at shell %g\n',TOL,shell_no);
  else
      disp('Direct sum did not converge!')
  end
end
fprintf('\n')
end

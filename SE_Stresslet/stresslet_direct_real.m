function [phi shellnorms] = stresslet_direct_real( idx, x, f, nvec, xi, L, nbox, TOL,varargin)
% Ewald summation for the stresslet -- Real space part.
%
% phi = stresslet_direct_real( idx, x, f, nvec, xi, L, nbox, TOL, [eval_x])
%        Evaluate potential, phi, at points x(idx) or eval_x. 
%        x    --  positions          N-by-3
%        nvec --  normal vectors   N-by-3
%        f    --  source strengths   N-by-3
%        xi   --  Ewald parameter
%        nbox --  periodic repications
%        TOL  -- desired tolerance
%        eval_x -- evaluate at points eval_x
%
% Example:
%   [x f] = generate_state(100,[1 1 1]);
%   stresslet_direct_real(1, x, f, nvec, 2, [1 1 1], 10)

VERBOSE = 0;
if ~exist('TOL','var')
    TOL=0;
end

p = -nbox:nbox; %boxes
[p1 p2 p3] = ndgrid(p,p,p);
p = [p1(:) p2(:) p3(:)];
Np = size(p,1);

pshell = max(abs(p'));
[pshell, I] = sort(pshell);
pshell = pshell';
p = p(I,:)*diag(L);
cprintf(VERBOSE,'\tComputing real space sum. Periodic images: %d\n', Np);
if numel(varargin)>0
    eval_x = varargin{1};
    noeval = size(eval_x,1);
    idx = 1:noeval;
else
    noeval = length(idx);  
    eval_x = x;
end
phi=zeros(noeval,3); 
shellnorms=zeros(nbox+1,noeval);
parfor ii=1:noeval
  m=idx(ii); 
  converged=0;
  shell_no=0;
  for jj=1:nbox+1    
    shell_no = jj-1;
    tmp = [0 0 0];
    shell_indices = find(pshell==shell_no);
    for j = shell_indices'
%     for j = 1:Np
      for n = 1:size(x,1) % particles
        if all(p(j,:)==0) && n==m % remove self interaction
          continue
        end
        r=eval_x(m,:)-x(n,:);
        tmp = tmp + (stresslet_op_real( r + p(j,:), nvec(n,:), xi)*f(n,:)')';
      end
    %    cprintf(mod(j,1000)==0,'\t%3d \t %5.2f%% \r',j, 100*j/Np);
    end
    phi(ii,:) = phi(ii,:) + tmp;
    shellnorm = norm(tmp);
    shellnorms(jj,ii)=shellnorm;
    cprintf(VERBOSE,'Shell %d: contrib=%g\n',shell_no,shellnorm);
    if shellnorm<TOL && shell_no>3
        converged=1;
        break
    end
  end
  if converged==1
      cprintf(VERBOSE,'Direct sum converged TOL=%g at shell %g\n',TOL,shell_no);
  else
      cprintf(VERBOSE,'Direct sum did not converge!\n')
  end
end
cprintf(VERBOSE,'\n')
end

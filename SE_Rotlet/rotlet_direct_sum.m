function u = rotlet_direct_sum( idx, x, t, L, nbox, TOL)

p = -nbox:nbox; %boxes
[p1 p2 p3] = ndgrid(p,p,p);
p = [L(1)*p1(:) L(2)*p2(:) L(3)*p3(:)];
Np = size(p,1);

pshell = max(abs(p'));
[pshell, I] = sort(pshell);
p = p(I,:);

fprintf('\tComputing direct sum. Periodic images: %d\n', Np);
noeval=length(idx);  
u=zeros(noeval,3); 
for ii=1:noeval
  m=idx(ii); 
  xm = x(m,:);
  converged=0;
  for shell_no=0:nbox    
    tmp = [0 0 0];
    parfor j = find(pshell==shell_no) % periodic images
        r = bsxfun(@plus, xm+p(j,:), -x);
        ri3 = sum(r.^2, 2).^(-3/2);
        txr = cross(t, r, 2);
        all_rotlets = bsxfun(@times, txr, ri3);
        if all(p(j,:)==0) % remove self interaction
            all_rotlets(m, :) = [];
        end
        tmp = tmp + sum(all_rotlets,1);
    end
    fprintf('Shell %d: contrib=%g\n',shell_no,norm(tmp));
    u(ii,:) = u(ii,:) + tmp;
    if norm(tmp)<TOL && shell_no>3
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

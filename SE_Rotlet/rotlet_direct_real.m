function [phi shellnorms] = rotlet_direct_real( idx, x, f, xi, L, varargin)
% Direct summation of rotlet real space part
%
% phi = rotlet_direct_real(idx, x, t, xi, box, 'layers', M, 'tol', TOL);
%   Compute interactions from all particles within M layers, stop at tolerance TOL.
%   
% phi = rotlet_direct_real(idx, x, t, xi, box, 'mode','cutoff', 'rc', rc);
%   Compute interactions from all particles within cutoff radius rc.

VERBOSE = 0;

p = inputParser();
p.addParameter('mode','layered');
p.addParameter('layers', 10);
p.addParameter('tol', 0);
p.addParameter('rc', inf);

p.parse(varargin{:});
nbox = p.Results.layers;
TOL = p.Results.tol;
rc = p.Results.rc;

switch p.Results.mode
  case 'layered'
    % default    
  case 'cutoff'
    assert(rc <= min(L), 'rc must be smaller than min(box)');
    nbox = 1;
    tol = 0;
  otherwise
    error('Unknown mode');
end
    
noeval = length(idx);  
eval_x = x;

p = -nbox:nbox; %boxes
[p1 p2 p3] = ndgrid(p,p,p);
p = [p1(:) p2(:) p3(:)];Np = size(p,1);

pshell = max(abs(p'));
[pshell, I] = sort(pshell);
pshell = pshell';
p = p(I,:)*diag(L);
cprintf(VERBOSE,'\tComputing real space sum. Periodic images: %d\n', Np);
phi=zeros(noeval,3); 
shellnorms=zeros(nbox+1,noeval);
rc2 = rc^2;
Npart = size(x,1);

perm = randperm(Npart);
x(perm,:) = x;
eval_x(perm,:) = eval_x;
f(perm,:) = f;
idx = perm(idx);

for ii=1:noeval
    m=idx(ii); 
    converged=0;
    shell_no=0;
    for jj=1:nbox+1    
        shell_no = jj-1;
        tmp = [0 0 0];
        shell_indices = find(pshell==shell_no);
        for j = shell_indices'            
            mask = true(Npart, 1); % interaction mask
            if all(p(j,:)==0) 
                mask(m) = false; % remove self interaction
            end
            r = bsxfun(@minus, eval_x(m,:) + p(j,:), x);
            r2 = sum(r.^2, 2);
            mask = mask & (r2 <= rc2);            
            if ~any(mask)
                continue
            end                       
            r2 = r2(mask);
            rn = sqrt(r2);
            fxr = cross(f(mask,:), r(mask,:), 2);
            A = ( erfc(xi*rn)./rn + 2*xi*exp(-xi^2*r2)/sqrt(pi) )./r2;            
            tmp = tmp + sum( bsxfun(@times, fxr, A), 1);
        end
        phi(ii,:) = phi(ii,:) + tmp;
        shellnorm = norm(tmp);
        shellnorms(jj,ii)=shellnorm;
        cprintf(VERBOSE,'Shell %d: contrib=%g\n',shell_no,shellnorm);
        if shellnorm<TOL && shell_no>1
            converged=1;
            break
        end
    end
    if isinf(rc)
        if converged==1
            cprintf(VERBOSE,'Real sum converged TOL=%g at shell %g\n',TOL,shell_no);
        else
            cprintf(VERBOSE,'Real sum did not converge!\n')
            %disp('Real sum did not converge!')
        end
    end
end


cprintf(VERBOSE,'\n')
end

function A = rotlet_op_real(x, xi)
    r2 = x(1)^2+x(2)^2+x(3)^2;
    r = sqrt(r2);
    c = xi^2*r2;
    A = x/r*(erfc(xi*r)/r2 + 2*xi*exp(-c)/(r*sqrt(pi)));
end
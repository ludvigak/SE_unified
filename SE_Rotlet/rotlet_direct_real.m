function [phi shellnorms] = rotlet_direct_real( idx, x, f, xi, L, varargin)
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
function [phi] = rotlet_real_rc(xe, x, f, xi, L, rc)

    assert(rc <= min(L), 'rc cannot exceed box size');

    nbox = 1;
    p = -1:1;
    [p1 p2 p3] = ndgrid(p,p,p);
    p = [p1(:)*L(1) p2(:)*L(2) p3(:)*L(3)];
    Np = size(p,1);

    Neval = size(xe, 1);
    Nsrc = size(x, 1);
    phi=zeros(Neval,3); 
    rc2 = rc^2;

    for ii=1:Neval
        r_base = bsxfun(@minus, xe(ii,:), x);
        tmp = [0 0 0];
        for j=1:Np
            mask = true(Nsrc, 1); % interaction mask
            r = bsxfun(@plus, r_base, p(j,:));
            r2 = sum(r.^2, 2);
            mask = mask & (r2 <= rc2) & (r2 > 0);
            if ~any(mask)
                continue
            end           
            r2 = r2(mask);
            rn = sqrt(r2);
            fxr = cross(f(mask,:), r(mask,:), 2);
            A = ( erfc(xi*rn)./rn + 2*xi*exp(-xi^2*r2)/sqrt(pi) )./r2;            
            tmp = tmp + sum( bsxfun(@times, fxr, A), 1);
        end
        phi(ii, :) = tmp;
    end
    
end

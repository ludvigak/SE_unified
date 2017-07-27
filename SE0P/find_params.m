function [M, P, rc, sl, nl, R] = find_params(L, Q, N, xi, tol)

[M, P, rc, ~] = ewald_k_inf2(L, Q, N, xi, tol*10);

for nl = 1:max(M/10,5);
    v = erfc(nl/2/xi); 
    if v<tol;
        break;
    end
end

for sl = 2:3
    if exp(-2*pi*nl*sl)<tol
        break
    end
end

R = L*sqrt(2)*(1+P/M/2);

for sl = 1:4
    v = 100/sqrt(M)*exp(-50/sqrt(M)*sl);
    if v<tol
        break
    end
end

for nl = 1:max(M/10,5);
    v = 2/sqrt(M)*exp(-160/M*nl); 
    if v<tol;
        break;
    end
end
sl = sl -1 ;


nl = -log10(tol);
if(tol>=1e-3)
    sl = 1;
elseif tol>=1e-6
    sl = 2;
elseif tol>=1e-12
    sl = 3;
else
    sl = 4;
end

nl = 1:20;
f = exp(-nl*M/(xi*L));
idx = find(f<tol);
nl = min(idx);

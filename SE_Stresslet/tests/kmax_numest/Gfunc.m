function [ G ] = Gfunc( j,l,m,kc,r,xi )

G = zeros(size(r));
for i=1:numel(r)
    G(i) = integral( @(k) kintegrand(j,l,m,k,r(i),xi), kc, inf, 'AbsTol',1e-3,'RelTol',1e-3);
end

end


function [k1 k2 k3] = k_vectors(M,box)
k1 = k_vec(M(1), box(1));
k2 = k_vec(M(2), box(2));
k3 = k_vec(M(3), box(3));

function [k] = k_vec(M,L)
if mod(M,2)==0
    MM = M/2;
    k = (2*pi/L)*[-MM:(MM-1)];
elseif mod(M-1,2)==0
    MM = (M-1)/2;
    k = (2*pi/L)*[-MM:MM];
else error('k-vectors not computed');
end

function [H info] = SE_grid(x, q, box, M, P, c)

% grid size
h = box(1)/M(1);

% input checking
assert(abs( h - box(2)/M(2)) < eps &&  abs( h - box(3)/M(3)) < eps);
assert(P <= min(M))
assert(mod(P-1,2) < eps)

H = zeros(M);
N = size(x,1);

for i=1:N
  [Q ax ay az] = SE_gaussian(x(i,:),h,M,P,c);
  
  % put into grid fcn
  H(ax,ay,az) = H(ax,ay,az) + q(i)*Q;
end

info.gaussian_mass_resid = abs(1-h^3*sum(Q(:)));
info.arithmetic_ratio = P^3/prod(M);
info.h = h;
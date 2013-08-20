function phi = SE_int(x, H, box, M, P, c)

% grid size
h = box(1)/M(1);

% input checking
assert(abs( h - box(2)/M(2)) < eps &&  abs( h - box(3)/M(3)) < eps);
assert(P <= min(M))
assert(mod(P-1,2) < eps)

[Q ax ay az] = SE_gaussian(x,h,M,P,c);

% integrate
phi = h^3*sum(sum(sum( Q.*H(ax,ay,az) )));

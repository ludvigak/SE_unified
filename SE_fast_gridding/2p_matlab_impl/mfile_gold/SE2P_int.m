function phi = SE2P_int(x, H, opt, c)

N = size(x,1);
phi = zeros(N,1);

for i = 1:N

  [Q ax ay az] = SE2P_gaussian(x(i,:),opt,c);

  phi(i) = opt.h^3*sum(sum(sum( Q.*H(ax,ay,az) )));

end

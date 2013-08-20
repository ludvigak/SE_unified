function [H info] = SE2P_grid(x, q, opt, c)

H = zeros([opt.M opt.M opt.Mz]);
N = size(x,1);

for i=1:N
    [Q ax ay az] = SE2P_gaussian(x(i,:),opt,c);
    
    % put into grid fcn
    H(ax,ay,az) = H(ax,ay,az) + q(i)*Q;
end

%info.gaussian_mass_resid = abs(1-opt.h^3*sum(Q(:)));
%info.arithmetic_ratio = opt.P^3/(opt.M^2*opt.Mz);
%info.h = opt.h;

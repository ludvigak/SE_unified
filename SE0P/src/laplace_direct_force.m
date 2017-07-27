function du = laplace_direct_force(x, f, box)

N = size(x, 1);

du = zeros(N, 3);
for target = 1:N
    source = [1:target-1 target+1:N];
    rvec   = bsxfun(@minus, x(target, :), x(source, :));
    dist   = sqrt(sum(rvec.^2, 2));
    u0     = f(source)./dist.^3;
    du(target,:) = sum(bsxfun(@times,u0,rvec));
end
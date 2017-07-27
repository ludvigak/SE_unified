function u = laplace_direct(x, f, box)

N = size(x, 1);

u = zeros(N, 1);
for target = 1:N
    source = [1:target-1 target+1:N];
    rvec = bsxfun(@minus, x(target, :), x(source, :));
    dist = sqrt(sum(rvec.^2, 2));
    u(target) = sum(f(source)./dist);
end

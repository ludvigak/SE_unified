clear
rng(1)
box = [2 2 2];
N = 100;

[x, f] = vector_system(N, box, 3);

udir = stokeslet_direct(x, f, box);

u = zeros(N, 3);
for i=1:N
    for j=[1:i-1 i+1:N]
        xij = x(i, :) - x(j, :);
        fj = f(j, :);
        r = norm(xij);
        u(i, :) = u(i, :) + ...
                  fj / r + xij * sum(xij.*fj) / r^3;
    end
end

error = udir - u;
max_err = norm(error(:), inf) / norm(u(:), inf)
assert(max_err < 1e-15, 'stokeslet direct failed')
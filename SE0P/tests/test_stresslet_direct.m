clear
rng(1)
box = [2 2 2];
N = 100;

[x, n] = vector_system(N, box, 3);
[~, q] = vector_system(N, box, 3);

udir = stresslet_direct(x, [q n], box);

u = zeros(N, 3);
for i=1:N
    for j=[1:i-1 i+1:N]
        xij = x(i, :) - x(j, :);
        nj = n(j, :);
        qj = q(j, :);        
        r = norm(xij);
        u(i, :) = u(i, :) - 6 * xij * sum(xij.*nj) * sum(xij.*qj) / r^5;
    end
end

error = udir - u;
max_err = norm(error(:), inf) / norm(u(:), inf);
assert(max_err < 1e-15, 'stresslet direct failed')

% Also test RS code with xi=0
opt.rc = inf;
opt.xi = 0;
urs = stresslet_real_space(x, [q n], opt);
error = udir - urs;
max_err = norm(error(:), inf) / norm(udir(:), inf);
assert(max_err < 2e-15, 'stresslet rs with xi=0 failed')

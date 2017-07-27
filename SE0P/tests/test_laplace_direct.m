clear
% Compare vectorized laplace_direct to naive sum

N = 100;
L = 2;
box = [L L L];
[x, f] = laplace_system(N, box);

u = laplace_direct(x, f, box);

uref = zeros(N, 1);
for i=1:N
    for j=[1:i-1 i+1:N]
        uref(i) = uref(i) + f(j)/norm(x(i,:)-x(j,:));
    end
end

error = u - uref;
rel_err_max = norm(error, inf) / norm(uref, inf);
assert(rel_err_max < 1e-14, 'laplace_direct failed')

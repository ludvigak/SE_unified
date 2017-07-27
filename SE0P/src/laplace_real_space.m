function u = laplace_real_space(x, f, opt)

xi = opt.xi;
rc = opt.rc;

N = size(x, 1);
[idx, d] = rangesearch(x, x, rc);

u = zeros(N, 1);
for target = 1:N
    nnb_idx = idx{target}(2:end);    
    nnb_r = d{target}(2:end); 
    u(target) = sum(f(nnb_idx)' .* erfc(xi*nnb_r)./nnb_r);
end

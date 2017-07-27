function du = laplace_real_space_force(x, f, opt)

xi = opt.xi;
rc = opt.rc;
c = 2/sqrt(pi)*xi;

N = size(x, 1);
[idx, d] = rangesearch(x, x, rc);


du = zeros(N, 3);
for target = 1:N
    nnb_idx= idx{target}(2:end);    
    nnb_r  = d{target}(2:end);
    r      = bsxfun(@minus,x(target,:),x(idx{target}(2:end),:));
    u1     = c*exp(-(xi^2*nnb_r.^2));
    u2     = erfc(xi*nnb_r)./nnb_r;
    u      = f(nnb_idx)'.*(u1+u2)./nnb_r.^2;
    du(target,:) = sum(bsxfun(@times,u',r));
end

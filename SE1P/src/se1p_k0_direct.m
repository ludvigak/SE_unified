function u = se1p_k0_direct(idx,x,f,opt)
    N = numel(f);
    xi = opt.xi;
    
    for target = idx        
        xt =  x(target,1:2);
        u_t = 0;
        for source = 1:N
            rho2 = norm(xt-x(source,1:2))^2;
            if(rho2>eps)
                v = rho2*xi*xi;
                u_t = u_t - f(source)*ein(v);
            end
        end
        u(target) = u_t/opt.box(3);
    end
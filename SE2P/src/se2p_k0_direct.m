function u = se2p_k0_direct(idx,x,f,opt,MATLAB)
    N = numel(f);
    xi = opt.xi;
    if(MATLAB)
        u = zeros(numel(idx),1);
        for target = idx        
            zt =  x(target,3);
            u_t = 0;
            for source = 1:N
                r = zt-x(source,3);
                u_t = u_t + f(source)*(exp(-xi^2*r^2)/xi+sqrt(pi)*r*erf(xi*r));
            end
            u(target) = -2*sqrt(pi)*u_t/opt.box(1)/opt.box(2);
        end
    else
        u = SE2P_direct_k0_mex(idx,x,f,opt);
    end
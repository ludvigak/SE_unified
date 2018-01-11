function u = se1p_k0_direct_force(idx,x,f,opt,varargin)
    N = numel(f);
    xi = opt.xi;
if nargin ==5
    MATLAB = varargin{1};
else
    MATLAB = false;
end    

if(MATLAB)
    nt = 1;
    u = zeros(numel(idx),3);
    for target = idx
        xt =  x(target,2:3);
        u_t = zeros(1,2);
        for source = 1:N
            dist = xt-x(source,2:3);
            rho2 = norm(dist)^2;
            if(rho2==0)
                continue;
            else
                v = rho2*xi*xi;
                u_t = u_t + f(source)*2.*dist/rho2*(1-exp(-v));
            end
        end
        
        u(nt,:) = -[0 u_t]/opt.box(1);
        nt = nt + 1;
    end
else
    u = SE1P_direct_k0_force_mex(idx,x,f,opt);
end
function u = se1p_fourier_space_direct(idx,x,f,opt) 

    N = numel(f);
    xi   = opt.xi;
    xi2        = xi*xi;
    TwoPiOverL = 2.*pi/opt.box(3);
    for target = idx        
        xt =  x(target,:);
        u_t = 0;
        for source = 1:N
            z = xt(3)-x(source,3);
            rho2 = norm(xt(1:2)-x(source,1:2))^2;
            b   = rho2*xi2;
            fn  = f(source);
            for p=1:opt.layers
                k3  = TwoPiOverL*p;
                a   = k3*k3/(4.*xi2);

                K0 = computeK0(a,b);

                k3z = -k3*z;
                u_t = u_t + 2*fn*cos(k3z)*K0;
            end
        end
        u(target) = u_t/(opt.box(3));
    end
end

function v = computeK0(a,b)
    func = @(t) exp(-a./t-b.*t)./t;
    v = integral(func,0,1,'reltol',1e-15);
end
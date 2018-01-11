function u = se1p_fourier_space_direct(idx,x,f,opt) 

    N = numel(f);
    xi   = opt.xi;
    xi2        = xi*xi;
    TwoPiOverL = 2.*pi/opt.box(1);
    for target = idx        
        xt =  x(target,:);
        u_t = 0;
        for source = 1:N
            X = xt(1)-x(source,1);
            rho2 = norm(xt(2:3)-x(source,2:3))^2;
            b   = rho2*xi2;
            fn  = f(source);
            for p=1:opt.layers
                k  = TwoPiOverL*p;
                a   = k*k/(4.*xi2);

                K0 = computeK0(a,b);

                kX = -k*X;
                u_t = u_t + 2*fn*cos(kX)*K0;
            end
        end
        u(target) = u_t/(opt.box(1));
    end
end

function v = computeK0(a,b)
    func = @(t) exp(-a./t-b.*t)./t;
    v = integral(func,0,1,'reltol',1e-15);
end
function [u varargout] = se1p_fourier_space_direct_force(idx,x,f,opt,varargin)
    
    if nargin == 5
        MATLAB = varargin{1};
    else
        MATLAB = false;
    end

    time = tic;
    if(MATLAB) 
        N = numel(f);
        u = zeros(numel(idx),3);
        xi   = opt.xi;
        xi2        = xi*xi;
        TwoPiOverL = 2.*pi/opt.box(3);
        nt = 1;
        for target = idx
            xt =  x(target,:);
            u_t = zeros(1,3);
            for source = 1:N
                rvec = xt(1:2)-x(source,1:2);
                rho2 = norm(rvec)^2;
                z = xt(3)-x(source,3);
                b   = rho2*xi2;
                fn  = f(source);
                for p=1:opt.layers
                    k3  = TwoPiOverL*p;
                    k3z = -k3*z;
                    a   = k3*k3/(4.*xi2);
                    
                    K0 = computeDbK0(a,b);
                    u_t(1:2) = u_t(1:2) - 4*fn*xi2*cos(k3z)*K0 * rvec(1:2);
                    
                    K0 = computeK0(a,b);
                    u_t(3)   = u_t(3) + 2*fn*k3*sin(k3z)*K0;
                end
            end
            u(nt,:) = u_t/(opt.box(3));
            nt = nt + 1;
        end
    else
        u = SE1P_direct_fd_force_mex(idx,x,f,opt);
    end
    walltime = toc(time);
    if(nargout==2)
        varargout{1} = walltime;
    end 

end

function v = computeK0(a,b)
    func = @(t) exp(-a./t-b.*t)./t;
    v = integral(func,0,1,'reltol',1e-15,'abstol',1e-15);
end

function v = computeDbK0(a,b)
    func = @(t) exp(-a./t-b.*t);
    v = integral(func,0,1,'reltol',1e-15,'abstol',1e-15);
end

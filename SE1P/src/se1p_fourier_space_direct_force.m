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
        TwoPiOverL = 2.*pi/opt.box(1);
        nt = 1;
        for target = idx
            xt =  x(target,:);
            u_t = zeros(1,3);
            for source = 1:N
                rvec = xt(2:3)-x(source,2:3);
                rho2 = norm(rvec)^2;
                X = xt(1)-x(source,1);
                b   = rho2*xi2;
                fn  = f(source);
                for p=1:opt.layers
                    k  = TwoPiOverL*p;
                    kX = -k*X;
                    a   = k*k/(4.*xi2);
                    
                    K0 = computeK0(a,b);
                    u_t(1)   = u_t(1) + 2*fn*k*sin(kX)*K0;
                    
                    K0 = computeDbK0(a,b);
                    u_t(2:3) = u_t(2:3) - 4*fn*xi2*cos(kX)*K0 * rvec;
                    
                end
            end
            u(nt,:) = u_t/(opt.box(1));
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

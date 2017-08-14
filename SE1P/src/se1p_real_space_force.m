function [u varargout]= se1p_real_space_force(idx, x, f, opt,varargin)

xi = opt.xi;

if nargin==6
    MATLAB = varargin{1};
    rc     = varargin{2};
elseif nargin ==5
    MATLAB = varargin{1};
    rc = inf;
else
    rc = inf;
    MATLAB = false;
end

c = 2/sqrt(pi)*xi;

N = size(x, 1);
time = tic;
if (MATLAB)
    u = zeros(numel(idx), 3);
    nt = 1;
    for target = idx
        for p=-opt.layers:opt.layers
            if p==0
                source = [1:target-1 target+1:N];
            else
                source = 1:N;
            end
            xsPv = x(source,:);
            xsPv(:,3) = xsPv(:,3) + p*opt.box(3);
            rvec = bsxfun(@minus, x(target, :), xsPv);
            dist = sqrt(sum(rvec.^2, 2));
            u1   = -c*exp(-(xi^2*dist.^2));
            u2   = -erfc(xi*dist)./dist;
            u0   = f(source).*(u1+u2)./dist.^2;

	    I = (dist<rc);
	    rvec = rvec(I,:);
	    u0 = u0(I);
	    if ~isempty(rvec)
                u(nt,:) = u(nt,:) + sum(bsxfun(@times,u0,rvec),1);
            end
        end
%        u(nt,:) = f(target) * u(nt,:);
        nt = nt + 1;
    end
else
    %   u = SE1P_direct_real_force_mex(idx,x,f,opt);
    u = SE1P_direct_rsrc_force_mex(idx,x,f,opt);
end
walltime = toc(time);

if(nargout==2)
   varargout{1} = walltime;
end 
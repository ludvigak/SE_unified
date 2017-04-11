function [u varargout]= se1p_direct_force(idx, x, f, opt, varargin)
 
N = size(x, 1);
if nargin==5
    MATLAB = varargin{1};
else
    MATLAB = true;
end

time = tic;
if(MATLAB)
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
            u0   = -f(source)./dist.^3;
            u(nt,:) = u(nt,:) + sum(bsxfun(@times,u0,rvec),1);
        end
%        u(nt,:) = f(target) * u(nt,:);
        nt = nt + 1;
    end
else
    u = SE1P_direct_force_mex(idx,x,f,opt);
end
walltime = toc(time);

if(nargout==2)
   varargout{1} = walltime;
end
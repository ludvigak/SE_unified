function [u varargout]= se1p_real_space(idx, x, f, opt,varargin)

xi = opt.xi;

if nargin==6
    MATLAB = varargin{1};
    rc     = varargin{2};
elseif nargin ==6
    MATLAB = varargin{1};
else
    rc = inf;
    MATLAB = false;
end

N = size(x, 1);
time = tic;
if (MATLAB)
    u = zeros(numel(idx), 1);
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
            if(dist<rc)
                u(target) = u(target) + sum(f(source) .* erfc(xi*dist)./dist);
            end
        end
    end
else
    u = SE1P_direct_real_mex(idx,x,f,opt);
end
walltime = toc(time);

if(nargout==2)
   varargout{1} = walltime;
end
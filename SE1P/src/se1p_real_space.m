function [u varargout]= se1p_real_space(idx, x, f, opt,varargin)

xi = opt.xi;

MATLAB = varargin{1};
rc = opt.rc;

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
            rvec = bsxfun(@minus, x(target,:),x(source,:));
            xsPv = rvec;
            xsPv(:,1) = xsPv(:,1) + p*opt.box(1);
            dist = sqrt(sum(xsPv.^2,2));

            I = (dist<rc);
            dist = dist(I);
            q = f(source(I));
            if ~isempty(dist)
                u(target) = u(target) + sum(q .* erfc(xi*dist)./dist);
            end
        end
    end
else
    u = SE1P_direct_rsrc_mex(idx,x,f,opt);
    %    u = SE1P_direct_real_mex(idx,x,f,opt);
end
walltime = toc(time);

if(nargout==2)
   varargout{1} = walltime;
end

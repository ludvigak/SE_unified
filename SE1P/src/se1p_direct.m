function [u varargout]= se1p_direct(idx, x, f, box, opt)
 
N = size(x, 1);
MATLAB = false;

time = tic;
if(MATLAB)
    u = zeros(N, 1);
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
            u(target) = u(target) + sum(f(source)./dist);
        end
    end
else
    u = SE1P_direct_mex(idx,x,f,opt);
end
walltime = toc(time);

if(nargout==2)
   varargout{1} = walltime;
end
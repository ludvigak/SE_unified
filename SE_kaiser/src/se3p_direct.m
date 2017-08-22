function [u varargout]= se3p_direct(idx, x, f,opt)

N = size(x, 1);

u = zeros(N, 1);
for target = idx
    for p1=-opt.layers:opt.layers
        for p2=-opt.layers:opt.layers
            for p3=-opt.layers:opt.layers
                if p1==0 && p2==0 && p3==0
                    source = [1:target-1 target+1:N];
                else
                    source = 1:N;
                end
                xsPv = x(source,:);
                xsPv = xsPv + [p1 p2 p3].*opt.box;
                rvec = bsxfun(@minus, x(target, :), xsPv);
                dist = sqrt(sum(rvec.^2, 2));
                u(target) = u(target) + sum(f(source)./dist);
            end
        end
    end
end

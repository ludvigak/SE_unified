function u = rotlet_direct(x, f, box)

N = size(x, 1);

MATLAB = 0;
if (MATLAB) 
    u = zeros(N, 3);
    for target = 1:N
        source = [1:target-1 target+1:N];
        rvec = bsxfun(@minus, x(target, :), x(source, :));
        ri3 = sum(rvec.^2, 2).^(-3/2);
        fxr = cross(f(source, :), rvec);
        u(target,:) = sum( bsxfun(@times, fxr, ri3) , 1);
    end
else
    u = SE0P_rotlet_direct_mex(x,f);
end

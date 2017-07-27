function u = stokeslet_direct(x, f, box)

N = size(x, 1);
MATLAB = 0;
if(MATLAB)
  u = zeros(N, 3);
  for target = 1:N
    source = [1:target-1 target+1:N];
    rvec = bsxfun(@minus, x(target, :), x(source, :));
    ri = 1./sqrt(sum(rvec.^2, 2));
    rdotf = sum(rvec .* f(source, :), 2);
    u(target,:) = sum( ...
		 bsxfun(@times, f(source, :), ri) + ...
		 bsxfun(@times, rvec, rdotf.*ri.^3) , 1);
  end
else
  u = SE0P_stokeslet_direct_mex(x,f);
end

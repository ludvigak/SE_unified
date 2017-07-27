function u = stresslet_direct(x, f, box)


q = f(:,[1 2 3]);
n = f(:,[4 5 6]);

N = size(x, 1);

MATLAB = false;

if(MATLAB)
  u = zeros(N, 3);
  for target = 1:N
    source = [1:target-1 target+1:N];
    rvec = bsxfun(@minus, x(target, :), x(source, :));
    ri5 = 1./sqrt(sum(rvec.^2, 2)).^5;
    rdotn = sum(rvec .* n(source, :), 2);
    rdotq = sum(rvec .* q(source, :), 2);
    u(target,:) = -6*sum( bsxfun(@times, rvec, rdotn.*rdotq.*ri5) , 1);
  end
else
  u = SE0P_stresslet_direct_mex(x,f);
end  

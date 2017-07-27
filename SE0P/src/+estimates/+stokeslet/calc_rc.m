function rc=calc_rc(f, opt, tol)
% Stokeslet real space cut-off
% rc = calc_rc(f, opt, tol)

Q = sum(f(:).^2);
V = prod(opt.box);
xi = opt.xi;
L = max(opt.box);


c = 2*sqrt(Q/V);

low = 0;
high = L;
for i=1:50
  rc = (low+high)/2;
  est = c*sqrt(rc)*exp(-(rc*xi)^2);
  if(est < tol)
    high = rc;
  else
    low = rc;
  end		       
end
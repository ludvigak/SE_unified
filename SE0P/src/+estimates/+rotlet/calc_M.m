function M = calc_M(f, opt, tol)
% rotlet grid size
% M = calc_M(f, opt, tol)

Q = sum(f(:).^2);
V = prod(opt.box);
xi = opt.xi;

c = sqrt(8*opt.xi^2*Q./(3*pi*V));
c = 1;
for M=2:2:100
  K = pi*M/opt.box(1);
  est =  c*1./sqrt(K).* exp(-K.^2/4/opt.xi^2);
  if(est<=tol)
     break;
  end
end
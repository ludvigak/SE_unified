function rc = calc_rc(f, opt,tol)
% stresslet real space cutoff
% rc = calc_rc(f, opt, tol)

q = f(:, [1 2 3]);
n = f(:, [4 5 6]);
S2 = zeros(3,3);
for l=1:3
    for m=1:3
        S2(l,m) = sum( ( q(:,l).*n(:,m) ).^2 );
    end
end   
Q = sum(S2(:));

V = prod(opt.box);
xi = opt.xi;

low = 0;
high = opt.box(1);
for i=1:50
  rc = (low+high)/2;
  xirc = xi*rc;
  % Hasimoto
  est = exp(-xirc.^2) * sqrt(112*xi^4*rc^3*Q/(9*V));
  % Beenaker
  % est = 4/3* exp(-xirc.^2) * sqrt(7*xi*xirc^3*Q/V);
  if(est < tol)
    high = rc;
  else
    low = rc;
  end		       
end






function P = calc_P(f, opt, tol)

xi = opt.xi;
Q = sum(f(:).^2);
c = sqrt(Q*xi*opt.box(1))/opt.box(1);

P = ceil(-2/(pi*.95^2)*log(tol/c));

if(mod(P,2))
    P = P +1;
end

P = min(P,34);
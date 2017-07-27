function est = real(f, opt)
% stresslet real space truncation error
% est = real(f, opt);

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
rc = opt.rc;

xirc = xi*rc;

% Hasimoto
% est = sqrt(Q/V*(...
%     -32/9*exp(-xirc^2)*sqrt(pi)*erfc(xirc) ...
%     +7*sqrt(2*pi)*erfc(sqrt(2)*xirc) ...
%     +16*pi*erfc(xirc)^2/(3*rc) ...
%     +4/27*exp(-2*xirc^2)*xi^2*rc*(-3+28*xirc^2)));
%est = 2*xi * exp(-xirc.^2) * sqrt(Q*rc*(28*xirc^2-3)/(9*V));
est = 4/3* exp(-xirc.^2) * sqrt(7*xi*xirc^3*Q/V);


est = exp(-xirc.^2) * sqrt(112*xi^4*rc^3*Q/(9*V));



% Beenakker
%est = sqrt( 3*exp(-2*xirc.^2) * Q/V * 1/27 * xi^2 .* rc .*( ...
%    327 + 1588*xirc.^2 -1392*xirc.^4 + 448*xirc.^6) );

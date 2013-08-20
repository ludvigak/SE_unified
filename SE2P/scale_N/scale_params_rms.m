function [xi M rc] = scale_params_rms(N,epsi)

% params
L = 1;
rho = N/(L^3);
g = N(1);
%C1 = 1e-3;
%C2 = 1e-2;
%epsi = 1e-10;

Q = N/4;

% verified formulas (wapr is a numerical approximation of lambertw)
%xxi = @(r) (1./(r*sqrt(2))) .* sqrt(wapr(2*C1^2./(r.^2*epsi^2)))
xxi = @(r) (1./r).*sqrt(wapr( (1/epsi)*sqrt(Q/(2*L^3)) ) );
kinf = @(x) sqrt(3)*L*x.*sqrt(wapr( 4*Q.^(2/3)./( 3*pi^(2/3)*L^2*x.^(2/3)*epsi^(4/3) ) ))/(2*pi);
%kinf = @(x) (x/pi)*sqrt(3/2) .* sqrt(wapr(2*pi^2*C2^(2/3)*x.^(2/3)./(3*epsi^(2/3))))

rc = (3*g./(4*pi*rho)).^(1/3);
xi = xxi(rc);
k = kinf(xi);
M = ceil(2*k);
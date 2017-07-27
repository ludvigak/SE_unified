function est = fourier_num(f, opt)

K = pi*opt.M(1)/opt.box(1);
Q = sum(f(:).^2);
V = prod(opt.box);
L = max(opt.box);
xi = opt.xi;
R = opt.R;

G = @(k) kernels.biharmonic(k, R);

B = @(k,r) G(k) .* 1./L^3 .* ...
    k.^2/4/xi^2 .* ...
    exp(-(k/2/xi).^2) .* ...
    2.*sin(k.*r)./(k.*r) ...
    .* k.^2 * sqrt(2)/3; % The directional average

% Integrate in k
errfunc = @(k,r) (L/2/pi)^3 * 2*pi*K.^2 .* B(k,r); 
delta = @(r) integral(@(k) errfunc(k,r), K, inf); 

% Integrate in r
Lmax = L/2;
est = sqrt(Q/L^3*integral(@(r) delta(r).^2.*r.^2*4*pi, 0, Lmax,'ArrayValued',true));

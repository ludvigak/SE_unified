function est = fourier(f, opt)
% stresslet fourier space truncation error
% est = fourier(f,opt)

q = f(:, [1 2 3]);
n = f(:, [4 5 6]);
S2 = zeros(3,3);
for l=1:3
    for m=1:3
        S2(l,m) = sum( ( q(:,l).*n(:,m) ).^2 );
    end
end   
Q = sum(S2(:));

K = pi*opt.M(1)/opt.box(1);
V = prod(opt.box);
xi = opt.xi;
L = max(opt.box);
R = opt.R;



est = sqrt(7*Q/6) * R * K^4 / (pi * xi^2 * L) * exp(-K^2/4/xi^2);

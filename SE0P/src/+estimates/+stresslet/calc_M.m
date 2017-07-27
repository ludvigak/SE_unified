function M = calc_M(f, opt, tol)
% stresslet grid size
% M = calc_M(f, opt, tol)

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
L = max(opt.box);

if(isfield(opt,'R'))
   R = opt.R;
else
   % R should be smaller than this value in 3d value since
   % R = |L_i+L_iP/2M_i|<|L_i+L_i/2|<=sqrt(3)*3*L/2;
   R = sqrt(3)*3*L/2;
end

c = sqrt(7*Q/6) * R / (pi * xi^2 * L);

M = 2:2:100;
K = M*pi/opt.box(1);
est = c*K.^4.*exp(-K.^2/(4*xi^2));
i = find(est>tol);
M = M(i(end)+1);
if(M>100)
  warning('FSE: GRID Calc Memory','The requested size results into a slow Fourier space sum; consider using a smaller Ewald parameter.');
end

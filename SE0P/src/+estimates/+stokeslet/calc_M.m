function M = calc_M(f, opt, tol)
% stokeslet grid size
% M = calc_M(f, opt, tol)

Q = sum(f(:).^2);
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

M = 2:2:140;
K = M*pi/opt.box(1);
est = sqrt(Q) * R * K.^3 / (pi * xi^2 * L) .* exp(-K.^2/4/xi^2);
i = find(est>tol)
M = M(i(end)+1);
if(M>100)
  warning('FSE: GRID Calc Memory','The requested size results into a slow Fourier space sum; consider using a smaller Ewald parameter.');
end


N = 100;

box = [1 2 3];
[x, t, xe] = generate_state(N, box);

xi = 8;
rc = 0.7;
ur = rotlet_real_rc(xe, x, t, xi, box, rc);
opt.P = 24;
opt.box = box;
opt.M = 40*box;
uk = SE_Rotlet(xe, x, t, xi, opt);
u1 = ur + uk;
opt.xi = xi;
opt.rc = rc;
u2 = rotlet_ewald_sum(xe, x, t, opt);

e = u1(:) - u2(:);
emax = norm(e(:), inf) / norm(u1(:), inf)
assert(emax < 1e-15);
disp('PASSED')


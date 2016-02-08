% Compare MEX and MATLAB codes for direct real space with cutoff radius

box = rand(1,3);

N = 50;
[x t xe] = generate_state(N, box);

opt.xi = 1;
opt.box = box;
opt.rc = min(box)/2;


f1 = @() rotlet_direct_rsrc(xe, x, t, opt);
f2 = @() rotlet_direct_real(1:size(xe,1), [xe;x], [xe*0; t], opt.xi, opt.box, 'mode', ...
                            'cutoff', 'rc', opt.rc);

u1 = f1();
u2 = f2();

e = u1 - u2;
emax = norm(e(:), inf) / norm(u1(:), inf);

assert(emax < 1e-15)
disp('PASSED')
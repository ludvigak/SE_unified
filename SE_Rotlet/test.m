%clear

rng(1);

N = 20/5;
x = rand(N,3);
t = rand(N,3);
%t(1,:) = 0

box = [1 1 1];
TOL = 1e-13;
idx = 1:N;
shells = 15;
k_shells = 15;


xi1 = 5;
xi2 = xi1+1;

%ud = rotlet_direct_sum(idx, x, t, box, shells, 1e-2)

ur1 = rotlet_direct_real(idx, x, t, xi1, box, 'layers',shells, 'tol', TOL);
uk1 = rotlet_direct_fd(idx, x, t, xi1, box, k_shells);
u1 = uk1+ur1

ur2 = rotlet_direct_real(idx, x, t, xi2, box, 'layers', shells, 'tol', TOL);
uk2 = rotlet_direct_fd(idx, x, t, xi2, box, k_shells);
u2 = uk2+ur2

%err=u2-ud
err=u2-u1

opt.P = 24;
opt.M = 31*box;
opt.box = box;

opt.M = 25*box;
ukSE1 = SE_Rotlet(idx, x, t, xi1, opt);
opt.M = 26*box;
ukSE2 = SE_Rotlet(idx, x, t, xi1, opt);

err1 = ukSE1-ukSE2
err2 = uk1-ukSE1

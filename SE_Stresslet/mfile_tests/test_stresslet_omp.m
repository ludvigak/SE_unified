function status = test_stresslet_omp()
rand('state',1)


N = 10000;
box = [1 2 3];
[x f n] = generate_state(N,box);
SE_opt.xi = 10;
SE_opt.box = box;
SE_opt.P = 16;
SE_opt.M = box*50;

SE_static  = SE_Stresslet_pre(x,SE_opt.xi,SE_opt);
u1 = SE_Stresslet(1:N, x, f, n, SE_opt.xi, SE_opt);
u2 = SE_Stresslet_omp(1:N, x, f, n, SE_opt.xi, SE_opt);
e = norm(u1(:)-u2(:),inf)/norm(u1(:),inf);

if e<1e-14
    fprintf('\n********** STRESSLET OMP: OK **********\n\n')
    status = 1;
else
    warning('STRESSLET OMP: FAILED')
    status = 0;
end
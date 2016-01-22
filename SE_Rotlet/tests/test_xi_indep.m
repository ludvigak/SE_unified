function status = test_xi_indep()

N = 100;

box = [1 2 3];
[x, t, xe] = generate_state(N, box);

u = {};
for xi = [8 9]
    rc = 0.7;
    ur = rotlet_real_rc(xe, x, t, xi, box, rc);
    opt.P = 24;
    opt.box = box;
    opt.M = 40*box;
    uk = SE_Rotlet(xe, x, t, xi, opt);
    u{end+1} = ur+uk;
end

err = u{1}-u{2};
max_err = norm(err(:), inf) / norm(u{1}(:), inf);

if max_err < 1e-14
    fprintf('\n********** XI INDEPENDENCE: OK **********\n\n')
    status = 1;
else
    warning('XI INDEPENDENCE: FAILED')
    status = 0;
end
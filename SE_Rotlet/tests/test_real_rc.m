function status = test_real_rc()

N = 100;

box = [1 2 3];
[x, t, xe] = generate_state(N, box);

idx = 1:N;
xe = x(idx,:);

xi = 8;
ref = rotlet_direct_real(idx, x, t, xi, box, 'layers', 1, 'tol', 1e-14);
rc = 0.75;
ur = rotlet_real_rc(xe, x, t, xi, box, rc);
max_err = norm(ref(:)-ur(:), inf) ./ norm(ref(:), inf);

if max_err < 1e-14
    fprintf('\n********** REAL RC: OK **********\n\n')
    status = 1;
else
    warning('REAL RC: FAILED')
    status = 0;
end
clear

N = 10;

box = [1 2 3];
[x, t, xe] = generate_state(N, box);

xi = 1;
k_shells =10;
ref = rotlet_direct_fd(xe, x, t, xi, box, k_shells);

results = [];
for P=16:32
    for M=10+(0:4)
        opt.P = P;
        opt.box = box;
        opt.M = M*box;
        uSE = SE_Rotlet(xe, x, t, xi, opt);
        err = norm(ref(:)-uSE(:), inf) / norm(ref(:), inf);
        results(end+1,:) = [P M log10(err)];
    end
end

if max(results(:,3)) < -10
    fprintf('\n********** FAST K-SPACE: OK **********\n\n')
else
    error('FAST K-SPACE: FAILED')
end
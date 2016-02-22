clear

N = 10;

box = [1 2 3];
[x, t, xe] = generate_state(N, box);

xi = 1;
results = [];
for P=16:32
    for M=10+(0:4)
        opt.P = P;
        opt.box = box;
        opt.M = M*box;
        opt.xi = xi;
        u = SE_Rotlet(xe, x, t, xi, opt);
        [U1, U2, U3] = SE_Rotlet(xe, x, t, xi, opt);
        F{1} = real( ifftn( ifftshift( U1 )));
        F{2} = real( ifftn( ifftshift( U2 )));
        F{3} = real( ifftn( ifftshift( U3 )));
        usplit = zeros(size(u));
        %usplit(:, 1) = SE_fg_int(xe, F1, opt);
        %usplit(:, 2) = SE_fg_int(xe, F2, opt);
        %usplit(:, 3) = SE_fg_int(xe, F3, opt);
        usplit = SE_fg_int(xe, F, opt);
        err = norm(u(:)-usplit(:), inf) / norm(u(:), inf);
        results(end+1,:) = [P M log10(err)];
    end
end

if max(results(:,3)) < -10
    fprintf('\n********** FAST K-SPACE: OK **********\n\n')
else
    error('FAST K-SPACE: FAILED')
end
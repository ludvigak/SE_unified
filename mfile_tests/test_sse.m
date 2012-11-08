function status = test_sse()

format short

% test SSE gridding / integration

box = [1 0.8 1.2]; % domain
N = 5000;        % numbver of charged particles
xi = 4;        % ewald parameter

% grid
SE_opt.M = 65*box
SE_opt.box = box;

% system
[x f nvec] = generate_state(N,box);
eval_idx = (1:N)';

% params
maxe = 0;

for P = [ 8 10 11 16 24];
    SE_opt.P = P;

    a=tic;
    SE_static  = SE_Stresslet_pre(x,xi,SE_opt);
    pretime = toc(a);

    [u_sse stats_sse]  = SE_Stresslet(eval_idx,x,f,nvec,xi,SE_opt,SE_static);
    [u stats]  = SE_Stresslet(eval_idx,x,f,nvec,xi,SE_opt);

    e = u-u_sse;
    e_inf = norm(e(:),inf)/norm(u(:),inf);

    pretime
    stats
    stats_sse

    fprintf('error = %.2e\n\n', e_inf);
    maxe = max(maxe, e_inf);
end
fprintf('max error = %.2e\n\n',maxe);
if maxe<1e-10
    fprintf('\n********** SSE CODE: OK **********\n\n')
    status = 1;
else
    warning('SSE CODE: FAILED')
    status = 0;
end
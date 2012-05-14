clear

% test SSE gridding / integration

box = [1 0.8 1.2]; % domain
N = 2e4;        % numbver of charged particles
xi = 4;        % ewald parameter

% grid
SE_opt.M = 50*box
SE_opt.box = box;

% system
[x f nvec] = generate_state(N,box);

% params
SE_opt.P = 16;
eval_idx = (1:N)';

a=tic;
SE_static  = SE_Stresslet_pre(x,xi,SE_opt);
pretime = toc(a);

[u_sse stats_sse]  = SE_Stresslet(eval_idx,x,f,nvec,xi,SE_opt,SE_static);
[u stats]  = SE_Stresslet(eval_idx,x,f,nvec,xi,SE_opt);

e = abs(u-u_sse)./abs(u);
disp('error=')
disp(e(1:10,:))

pretime
stats
stats_sse

maxe = max(e(:))

if maxe<1e-10
    fprintf('\n********** SSE CODE: OK **********\n\n')
else
    error('SSE CODE: FAILED')
end
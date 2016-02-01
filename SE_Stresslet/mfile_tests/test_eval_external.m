clear

rand('state',1)

box = [1 0.5 0.1]; % domain
N = 100;        % numbver of charged particles
xi = 4;        % ewald parameter

% grid
SE_opt.M = 30*box
SE_opt.box = box;

% system
[x f nvec] = generate_state(N,box);

% params
SE_opt.P = 24;
eval_idx = 1:N;
kmax = floor( (SE_opt.M(1)-1)/2 );

u  = SE_Stresslet(eval_idx,x,f,nvec,xi,SE_opt);

SE_opt.eval_external = 1;
SE_opt.eval_x = x(eval_idx, :);
u_ext  = SE_Stresslet(eval_idx,x,f,nvec,xi,SE_opt);

res = max(abs(u - u_ext))

if res<1e-10
    fprintf('\n********** EVAL EXTERNAL: OK **********\n\n')
else
    error('EVAL EXTERNAL: FAILED')
end


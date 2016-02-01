clear
rand('state',1)

box = [1 0.5 0.1]; % domain
N = 10;        % numbver of charged particles
xi = 4;        % ewald parameter

% grid
SE_opt.M = 30*box
SE_opt.box = box;

% system
[x f nvec] = generate_state(N,box);

% params
SE_opt.P = 24;
eval_idx = 1;
kmax = floor( (SE_opt.M(1)-1)/2 );

u_fast  = SE_Stresslet(eval_idx,x,f,nvec,xi,SE_opt)
% 
u_dir = stresslet_direct_fd( eval_idx, x, f, nvec, xi, SE_opt.box, kmax)

res = max(abs(u_fast-u_dir))

if res<1e-10
    fprintf('\n********** FAST SUMMATION: OK **********\n\n')
else
    error('FAST SUMMATION: FAILED')
end
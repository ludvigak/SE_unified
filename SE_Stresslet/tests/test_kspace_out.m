clear

format short

box = [1 0.8 1.2]; % domain
N = 5000;        % numbver of charged particles
xi = 4;        % ewald parameter

% grid
SE_opt.M = 60*box
SE_opt.box = box;
SE_opt.xi = xi;

% system
[x f nvec] = generate_state(N,box);
eval_idx = (1:N)';

% params
maxe = 0;

for P = [ 25 8 12 16 24];
    SE_opt.P = P;

    u  = SE_Stresslet(eval_idx,x,f,nvec,xi,SE_opt);
    [G1,G2,G3]  = SE_Stresslet(eval_idx,x,f,nvec,xi,SE_opt);
    
    F{1} = real( ifftn( ifftshift( G1 )));
    F{2} = real( ifftn( ifftshift( G2 )));
    F{3} = real( ifftn( ifftshift( G3 )));
    uf = SE_fg_int(x, F, SE_opt);
    
    e = u-uf;
    e_inf = norm(e(:),inf)/norm(u(:),inf);
    assert(e_inf < 2e-14, 'FOURIER OUT: FAILED');
end

fprintf('\n********** FOURIER OUT: OK **********\n\n')

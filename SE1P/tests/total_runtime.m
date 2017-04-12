clear
rng(1);

N = 20000;
L = 5;
box = [L L L];
[x, f] = vector_system(N, box);

opt.M = 46;
opt.xi = .1;
opt.P = 16;
opt.rc = 2.5;
opt.box = box;
opt.layers = 0;
opt.sl = 3;
opt.nl = 4;
opt.s0 = 2.5;


%% Ewald
disp('One periodic Ewald...')
[uf tfourier]= se1p_fourier_space_force(x, f, opt);
[ur treal]= se1p_real_space(1:N,x, f, opt,false);
tic
ur3 = SE1P_rsrc_cell_mex(x', f, opt.rc, opt.xi, opt.layers, box);
toc
treal
norm(ur3-ur)/N

% [ur treal]= se1p_real_space_force(1:N,x, f, opt,false);
% ur3 = SE1P_rsrc_cell_force_mex(x', f, opt.rc, opt.xi, opt.layers,box)';
% norm(ur3-ur)/N
return
ue = uf+ur;

%% Direct
% idx = 1:N;
% disp('Direct sum...')
% opt.P=24;
% u1 = se1p_fourier_space_direct_force(idx,x, f, opt, false);
% u2 = se1p_real_space_force(idx,x, f, opt,true);
% u3 = se1p_k0_direct_force(idx,x,f,opt,true);
% u  = u1+u2+u3;

% ue = ue(idx,:);
% rms_error = rms(u(:)-ue(:)) / rms(u(:))

time = tfourier.total+treal
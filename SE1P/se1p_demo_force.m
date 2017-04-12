clear
rng(1);

N = 20;
L = 1;
box = [L L L];
[x, f] = vector_system(N, box);

M0 = 32; % Set M0 to an even number, the rest is auto

opt.M = M0*box(1);
opt.xi = pi*M0 / 12;
opt.P = 28;
opt.rc = 6 / opt.xi; opt.rc = min(opt.rc,L);
opt.box = box;
opt.layers = 15;
opt.sl = 3;
opt.nl = 4;
opt.s0 = 2.5;


%% Ewald
disp('One periodic Ewald...')
[uf tfourier]= se1p_fourier_space_force(x, f, opt);
[ur treal]= se1p_real_space_force(1:N,x, f, opt,true);
uk = se1p_k0_direct_force(1:N,x,f,opt,true);
ue = uf+ur;

%% Direct
disp('Direct sum...')
[u tdirect]= se1p_direct_force(1:N, x, f, opt,true);
opt.P=32;
[u1 tfourier]= se1p_fourier_space_direct_force(1:N,x, f, opt, false);
[u2 treal]= se1p_real_space_force(1:N,x, f, opt,true);
u3 = se1p_k0_direct_force(1:N,x,f,opt,true);
u = u1+u2+u3;
rms_error = rms(u(:)-ue(:)) / rms(u(:))

% fprintf('fourier-sp time\t real-sp time\t total time\t direct time\n');
% fprintf('%10.6f\t%10.6f\t%10.6f\t%10.6f\n',...
%        tfourier.total,treal,tfourier.total+treal,tdirect);
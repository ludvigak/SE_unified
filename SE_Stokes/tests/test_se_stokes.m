clear all

box = [1 1 1]; % domain
N = 10;        % numbver of charged particles
xi = 4;        % ewald parameter

SE_opt.P = 24; % support 

% grid
SE_opt.M = 31*box;
SE_opt.box = box;

% system
[xs fs] = SE_charged_system(N,box,'vector');
[xe fe] = SE_charged_system(N,box,'vector');
fe(:) = 0;

x = [xs; xe];
f = [fs; fe];

% parameters for (reference) direct Ewald sum
ED_opt.layers = (SE_opt.M(1)-1)/2;
ED_opt.xi = xi;
ED_opt.box = box;

% compute FD Ewald sum
ref = SE3P_Stokes_direct_fd_mex(1:2*N,x,f,ED_opt);
refnorm = norm(ref(:), inf);

u1 = SE_Stokes(1:2*N,x,f,xi,SE_opt);
e1 = norm(u1(:) - ref(:), inf) / refnorm;
assert(e1 < 1e-13, 'SE_Stokes differed from direct code')

u2 = SE_Stokes_par(1:2*N,x,f,xi,SE_opt);
e2 = norm(u2(:) - ref(:), inf) / refnorm;
assert(e2 < 1e-13, 'SE_Stokes_par differed from direct code')

u3 = SE_Stokes_ext(xe,xs,f,xi,SE_opt);
e3 = norm(reshape(u3 - ref(N+1:end,:), [], 1), inf) / refnorm;
assert(e3 < 1e-13, 'SE_Stokes_ext differed from direct code')

disp('PASSED')
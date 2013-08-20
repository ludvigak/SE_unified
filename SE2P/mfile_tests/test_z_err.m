clear all,  close all

SE_opt.box = [1 1 1]; % domain
N = 200;              % number of charged particles
xi = 8;               % ewald parameter

SE_opt.P = 25;
s = [2 3 4];

% grid
SE_opt.M = 35;

% charge-neutral system
[x q] = SE_charged_system(N,SE_opt.box,'scalar');

% parameters for (reference) direct Ewald sum
ED_opt.layers = (SE_opt.M-1)/2;
ED_opt.xi = xi;
ED_opt.box = SE_opt.box;

% compute FD Ewald sum
ref = SE2P_direct_fd_mex(1:N,x,q,ED_opt);

% fast method
sty = {'.','*','+'};
hold on
for i = 1:length(s)
    % for oversampling factor, s
    SE_opt.sampling_factor=s(i);

    u = spectral_ewald_2P(1:N,x,q,xi,SE_opt);
    e = abs(u - ref);
    plot(x(:,3),e,sty{i})
end
publication_fig
set(gca,'YScale','log')
xlabel('z position'), ylabel('abs error')
set(gca,'YTick',[1e-15 1e-10 1e-5 1e-0])
axis([0 1 1e-16 1e0])
fname = sprintf('output/SE_z_accuracy_xi%d_M%d_P%d',xi,SE_opt.M(1),SE_opt.P);
write_fig(1,fname)
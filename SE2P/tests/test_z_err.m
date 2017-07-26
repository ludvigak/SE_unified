clear all,  close all

opt.box = [1 1 1]; % domain
N = 200;              % number of charged particles
xi = 8;               % ewald parameter

opt.P = 25;
s = [2 3 4];
opt.s0=3;

% grid
opt.M = 35;

% charge-neutral system
[x q] = vector_system(N,opt.box);

% parameters for (reference) direct Ewald sum
opt.layers = (opt.M-1)/2;
opt.xi = xi;

% compute FD Ewald sum
ref = SE2P_direct_fd_mex(1:N,x,q,opt);
ref0= SE2P_direct_k0_mex(1:N,x,q,opt);
ref = ref+ref0;

% fast method
sty = {'.','*','+'};
hold on
for i = 1:length(s)
    % for oversampling factor, s
    opt.s=s(i);

    u = se2p_fourier_space(x,q,opt);
    e = abs(u - ref);
    plot(x(:,3),e,sty{i})
end
publication_fig
set(gca,'YScale','log')
xlabel('z position'), ylabel('abs error')
set(gca,'YTick',[1e-15 1e-10 1e-5 1e-0])
axis([0 1 1e-16 1e0])
%fname = sprintf('output/SE_z_accuracy_xi%d_M%d_P%d',xi,opt.M(1),opt.P);
%write_fig(1,fname)

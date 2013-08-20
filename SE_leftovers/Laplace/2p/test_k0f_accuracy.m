clear, close all

N = 10;
box = [1 1 1];
%opt.layers=1; % FIXME: this prevents a core dump!
opt.xi=5;
opt.box = box;
opt.method='mex';
%opt.k0_M=40;
opt.interp_meth = 'spectral';
[x q] = SE_state(N,box,1);

tic
phi_1 = ewald_direct_2dp_k0_mex(1:N, x, q, opt);
toc

M = 2:2:63;
for j=1:length(M)
    opt.k0_M=M(j);
    opt.interp_meth = 'spectral';
    %tic
    phi_2 = ewald_direct_2dp_k0_fast(1:N, x, q, opt);
    %toc

    einf_1(j) = max(abs(phi_1-phi_2));
    
    %opt.interp_meth = 'spline';
    %phi_2 = ewald_direct_2dp_k0_fast(1:N, x, q, opt);
    %einf_2(j) = max(abs(phi_1-phi_2));
end

xi = opt.xi;
est = xi^3*2.^(-pi*xi^(-1/2)*(M-1));

semilogy(M,einf_1,'.-', M,est,'-r')

publication_fig, set(gca,'YLim',[1e-16 1]), grid on
xlabel('P')

fname = sprintf('output/convg_k0_cheb_xi%d_N%d',opt.xi,N)
write_fig(1,fname)
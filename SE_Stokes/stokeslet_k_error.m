clear
rng(1);

L = 2;
xi = 10;
N = 2000;
box = L*[1 1 1];
[x f] = SE_charged_system(N,box,'vector');
M_max = round(14*xi/pi)
M_ref = round(15*xi/pi);
opt.P = 32;
opt.box = box;
opt.M = M_ref*box;
ref = SE_Stokes(1:N, x, f, xi, opt);
ref_max = norm(ref(:), inf);

Mlist = 2:2:M_max;
[err_inf, err_rms, est, est_L, est_num] = deal(zeros(size(Mlist)));
tic
for i=1:numel(Mlist)
    this_opt = opt;
    this_opt.M = Mlist(i)*box;    
    uk = SE_Stokes(1:N, x, f, xi, this_opt);
    err = uk - ref;
    err_rms(i) = sqrt(1/N*sum(err(:).^2));
    est(i) =  stokeslet_est_k(f, this_opt, xi);
    %est_L(i) =  stokeslet_est_k_Lindbo(f, this_opt, xi);
    %est_num(i) =  stokeslet_est_k_numint(f, this_opt, xi);
end
toc
K = Mlist*pi;
clf
semilogy(K / xi, err_rms / ref_max, '.', 'DisplayName','Measured')    
hold on
semilogy(K / xi, est / ref_max,'-', 'DisplayName','Estimate')
%semilogy(K / xi, est_L / ref_max,'-', 'DisplayName','Lindbo')
%semilogy(K / xi, est_num / ref_max,'-', 'DisplayName','Numerical')

ylim([1e-16 1])
grid on
xlabel('$K/\xi$','interpreter','latex')
ylabel('Stokeslet Fourier space RMS error (rel.)')
title(sprintf('N=%d, \\xi=%g, L=%d',N,xi,L)) 
legend toggle
drawnow

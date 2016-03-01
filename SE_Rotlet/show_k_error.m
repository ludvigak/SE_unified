clear

L = 1;
xi = 20;
N = 1000;
box = L*[1 1 1];
[x t xe] = generate_state(N, box);
M_max = round(14*xi/pi);
M_ref = round(15*xi/pi);
opt.P = 32;
opt.box = box;
opt.M = M_ref*box;
ref = SE_Rotlet(xe, x, t, xi, opt);
ref_max = norm(ref(:), inf);

Mlist = 2:2:M_max;
err_inf = [];
err_rms = [];
est = [];
for i=1:numel(Mlist)
    this_opt = opt;
    this_opt.M = Mlist(i)*box;    
    uk = SE_Rotlet(xe, x, t, xi, this_opt);
    err = uk - ref;
    err_rms(i) = sqrt(1/N*sum(err(:).^2));
    est(i) =  rotlet_est_k(t, this_opt, xi);
end

K = Mlist*pi;
clf
semilogy(K / xi, err_rms / ref_max, '.')    
hold on
semilogy(K / xi, est / ref_max,'-')
ylim([1e-16 1])
grid on
xlabel('$K/\xi$','interpreter','latex')
ylabel('Fourier space RMS error (rel.)')
legend('Measured','Estimate');
drawnow

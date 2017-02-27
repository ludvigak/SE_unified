clear

L = 1;
xi = 20;
N = 1000;
ref_shells = 1;
box = L*[1 1 1];
rc_max = 7/xi;
rc_list = linspace(0,rc_max, 50);
rc_list = rc_list(2:end);
[x t] = generate_state(N, box);
idx = 1:N;
opt.box = box;
opt.xi = xi;

tic
%ur_ref = rotlet_direct_real(idx, x, t, xi, box, 'layers', ref_shells, 'tol', 1e-14);
opt.rc = min(box);
p = randperm(N); % Permute to actually get any round-off errors
ur_ref = rotlet_direct_rsrc(x(idx,:), x(p,:), t(p,:), opt);
toc
ref_max = norm(ur_ref(:), inf);

start = tic();
for i=1:numel(rc_list);
    opt.rc = rc_list(i);
    ur = rotlet_direct_rsrc(x(idx,:), x, t, opt);
    err = abs(ur - ur_ref);
    err_rms(i) = sqrt(1/N*sum(err(:).^2)) / ref_max;
end
toc(start);

clf
semilogy(xi*rc_list,err_rms,'.')
hold on
rc = rc_list;
semilogy(xi*rc, rotlet_est_real(t, box, rc, xi) / ref_max);
semilogy(xi*rc, rotlet_est_real_long(t, box, rc, xi) / ref_max, '--')
ylim([min(err_rms(err_rms>0))/2 1]);
xlim([0 max(xi*rc)]);
grid on
ylabel('Real space RMS error (rel.)')
xlabel('$\xi r_c$','interpreter','latex');
legend('Measured','Estimate (simple)','Estimate (full)')

clear

L = 2;
xi = 5;
N = 200;
ref_shells = 1;
box = L*[1 1 1];
rc_list = linspace(0,min(box)*0.75, 50);
rc_list = rc_list(2:end);
[x t] = generate_state(N, box);
idx = 1:N;

tic
ur_ref = rotlet_direct_real(idx, x, t, xi, box, 'layers', ref_shells, 'tol', 1e-14);
toc
ref_max = norm(ur_ref(:), inf);

start = tic();
parfor i=1:numel(rc_list);
    rc = rc_list(i);
    ur = rotlet_real_rc(x(idx,:), x, t, xi, box, rc);
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
ylabel('Real space RMS error (rel.)')
xlabel('$\xi r_c$','interpreter','latex');
legend('Measured','Estimate (simple)','Estimate (full)')

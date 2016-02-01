clear

L = 2;
xi = 8;
N = 2000;
ref_shells = 1;
box = L*[1 1 1];
rc_list = linspace(0,min(box)*0.49, 100);
rc_list = rc_list(2:end);
[x f n] = generate_state(N, box);
idx = 1:N;

tic
ur_ref = stresslet_direct_real_fast(idx, x, f, n, xi,  box, ref_shells, L);
toc
ref_max = norm(ur_ref(:), inf);

start = tic();
parfor i=1:numel(rc_list);
    rc = rc_list(i);
    ur = stresslet_real_rc( x, f, n, xi, box, rc);
    err = abs(ur - ur_ref);
    err_rms(i) = sqrt(1/(3*N)*sum(err(:).^2)) / ref_max;
end
toc(start);

S2 = zeros(3,3);
for l=1:3
    for m=1:3
        S2(l,m) = sum( ( f(:,l).*n(:,m) ).^2 );
    end
end   
Q = sum(S2(:));
V = prod(box);
rc = rc_list;
xirc = xi*rc;
est = sqrt( exp(-2*xirc.^2) * Q/V * 1/27 * xi^2 .* rc .*( ...
    327 + 1588*xirc.^2 -1392*xirc.^4 + 448*xirc.^6) );

clf
semilogy(xi*rc_list,err_rms, '.')
hold on
semilogy(xi*rc, est / ref_max);
%semilogy(xi*rc, rotlet_est_real_long(t, box, rc, xi) / ref_max, '--')
ylim([min(err_rms(err_rms>0))/2 1]);
xlim([0 max(xi*rc)]);
grid on
ylabel('Real space RMS error (rel.)')
xlabel('$\xi r_c$','interpreter','latex');
legend('Measured','Estimate')

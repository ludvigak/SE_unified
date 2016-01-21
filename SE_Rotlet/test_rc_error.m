clear


L = 2;
xi = 5;
N = 300;
ref_shells = 1;

box = L*[1 1 1];
x = L*rand(N,3);
t = 1-2*rand(N,3);
idx = 1:N;



ur_ref = rotlet_direct_real(idx, x, t, xi, box, 'layers', ref_shells, 'tol', 1e-14);
ref_max = norm(ur_ref(:), inf);

rc_list = linspace(0,min(box)*0.75, 20);


start = tic();
parfor i=1:numel(rc_list);
    rc = rc_list(i);
    ur = rotlet_direct_real(idx, x, t, xi, box, 'mode', 'cutoff', 'rc', rc);
    err = abs(ur - ur_ref);
    err_rms(i) = sqrt(1/N*sum(err(:).^2)) / ref_max;
    err_inf(i) = norm(err(:),inf) / ref_max;
end
toc(start);

clf
semilogy(rc_list,err_rms,'o-')
grid on
hold on
%semilogy(rc_list,err_inf,'d-')

rc = rc_list;
Q = sum(t(:).^2);
V = prod(box);

%semilogy(rc, sqrt(2/3)*2*sqrt(Q/V./rc).*exp(-(xi*rc).^2) / ref_max);
semilogy(rc, sqrt(8*Q./(3*V*rc)).*exp(-(xi*rc).^2) / ref_max);

%est2 = erfc(xi*rc).^2./rc + sqrt(2/pi)*xi*erfc(sqrt(2)*rc*xi);
%est = sqrt(Q/V*4*pi*est2);
%semilogy(rc, est / ref_max, '--')

ylim([min(err_rms(err_rms>0))/2 1]);
xlim([0 max(rc)]);
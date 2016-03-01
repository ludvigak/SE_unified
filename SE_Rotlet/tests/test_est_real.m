clear

L = 1;
xi = 6;
N = 400;
ref_shells = 1;
box = L*[1 1 1];

xirc = linspace(1,5,10);
rc_list = xirc/xi;
[x t] = generate_state(N, box);
xe = x;

opt.xi = xi;
opt.box = box;
opt.rc = min(box);

p = randperm(N); % Permute to actually get any round-off errors
ur_ref = rotlet_direct_rsrc(xe, x(p,:), t(p,:), opt);
ref_max = norm(ur_ref(:), inf);


for i=1:numel(rc_list);
    opt.rc = rc_list(i);
    ur = rotlet_direct_rsrc(xe, x, t, opt);
    err = abs(ur - ur_ref);
    err_rms(i) = sqrt(1/N*sum(err(:).^2));
end

rc = rc_list;

est_short = rotlet_est_real(t, box, rc, xi);
est_long = rotlet_est_real_long(t, box, rc, xi);

err_short = max(abs(1-abs(err_rms ./ est_short)));
err_long = max(abs(1-abs(err_rms ./ est_long)));

short_passed = err_short < 0.3;
long_passed = err_long < 0.15;

if short_passed && long_passed
    fprintf('\n********** EST REAL: OK **********\n\n')
else
    error('EST REAL: FAILED')
end    

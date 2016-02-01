clear

L = 1;
xi = 15;
N = 400;
box = L*[1 1 1];
[x t xe] = generate_state(N, box);
M_min = ceil(4*xi/pi);
M_max = floor(10*xi/pi);
M_ref = round(12*xi/pi);
opt.P = 32;
opt.box = box;
opt.M = M_ref*box;
ref = SE_Rotlet(xe, x, t, xi, opt);
ref_max = norm(ref(:), inf);

Mlist = M_min:1:M_max;
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

q=abs(1-abs(est./err_rms));
if q < .5
    fprintf('\n********** EST K: OK **********\n\n')
else
    error('EST K: FAILED')
end    

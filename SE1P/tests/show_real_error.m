function output = show_real_error(varargin)

p = inputParser;
p.addParameter('L',2);
p.addParameter('N',100);

p.parse(varargin{:});
input = p.Results;
L = input.L;
N = input.N;


rng(1);
box = L*[1 1 1];
[x f] = vector_system(N, box);f = (-1).^(1:N)';
rcmax = L/2;
opt.xi = 7 / rcmax;
opt.box = box;
opt.layers = 100;

% compute reference
opt.rc = inf;
perm = randperm(N); % permute to show roundoff errors
ref = SE1P_direct_rsrc_mex(1:N, x(perm,:), f(perm,:), opt);
ref(perm,:) = ref;
ref_max = rms(ref);
% Iterate
rclist = linspace(0,rcmax, 30);rclist = rclist(2:end);
[err_inf, err_rms, est] = deal(zeros(size(rclist)));
tic
for i=1:numel(rclist)
    this_opt = opt;
    this_opt.rc = rclist(i);
    ur = SE1P_direct_rsrc_mex(1:N, x, f, this_opt);
    err = ur - ref + eps(ref_max);
    err_rms(i) = rms(err);
    est(i) = real_estimate(f, this_opt);
end
toc

ref_max=1;
clf
xirc = opt.xi * rclist;
semilogy(xirc, err_rms / ref_max, 'b.', 'DisplayName','Measured')    
hold on
semilogy(xirc, est / ref_max,'r-', 'DisplayName','Estimate')
ylim([1e-17 1])
grid on
xlabel('$\xi r_c$','interpreter','latex')
ylabel('real space RMS error (rel.)')
title(sprintf('N=%d, \\xi=%g, L=%d',N,opt.xi,L)) 
legend toggle
drawnow


output.xirc = xirc;
output.err_rms_rel = err_rms; %/ ref_max;
output.est_rel = est;% / ref_max;
output.xi = opt.xi;


function e = real_estimate(f, opt)
Q = sum(f(:).^2);
xi = opt.xi;
L = opt.box(1);
rc = opt.rc;

e = sqrt(Q*rc/(2*L^3))*(xi*rc)^(-2)*exp(-xi^2*rc^2);

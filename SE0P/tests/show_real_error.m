function output = show_real_error(kernel, varargin)
% show_real_error(kernel, [args])
%
% ex. show_fourier_error('stokeslet')

p = inputParser;
p.addParameter('L',2);
p.addParameter('N',100);

p.parse(varargin{:});
input = p.Results;
L = input.L;
N = input.N;


fse = fsewald(kernel);

rng(1);
box = L*[1 1 1];
[x f] = fse.generator(N, box);
x(1,:) = 0; x(2,:) = box*(1-eps);
rcmax = L/2;
opt.xi = 7 / rcmax;
opt.box = box;

% compute reference
opt.rc = inf;
perm = randperm(N); % permute to show roundoff errors
ref = fse.real_sum(x(perm,:), f(perm,:), opt);
ref(perm,:) = ref;
ref_max = norm(ref(:), inf);
% Iterate
rclist = linspace(0,rcmax, 50);rclist = rclist(2:end);
[err_inf, err_rms, est, err_3p_rms] = deal(zeros(size(rclist)));
tic
for i=1:numel(rclist)
    this_opt = opt;
    this_opt.rc = rclist(i);
    ur = fse.real_sum(x, f, this_opt);    
    err = ur - ref + eps(ref_max);
    err_rms(i) = sqrt(1/N*sum(err(:).^2));
    est(i) = fse.real_estimate(f, this_opt);
end
toc

clf
xirc = opt.xi * rclist;
semilogy(xirc, err_rms / ref_max, '.', 'DisplayName','Measured')    
hold on
semilogy(xirc, est / ref_max,'-', 'DisplayName','Estimate')
ylim([1e-17 1])
grid on
xlabel('$\xi r_c$','interpreter','latex')
ylabel(sprintf('%s real space RMS error (rel.)', kernel))
title(sprintf('N=%d, \\xi=%g, L=%d',N,opt.xi,L)) 
legend toggle
drawnow


output.xirc = xirc;
output.err_rms_rel = err_rms / ref_max;
output.est_rel = est / ref_max;
output.xi = opt.xi;
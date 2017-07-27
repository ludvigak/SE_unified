function output = show_fourier_error(kernel, varargin)
% show_fourier(kernel, [args])
%
% ex. show_fourier_error('stokeslet')


p = inputParser;
p.addParameter('M',20);
p.addParameter('L',2);
p.addParameter('N',100);
p.addParameter('do_3p',false);

p.parse(varargin{:});
input = p.Results;
M_max = input.M;
do_3p = input.do_3p;
L = input.L;
N = input.N;

fse = fsewald(kernel);

box = L*[1 1 1];
[x f] = fse.generator(N, box);
x(1,:) = 0; x(2,:) = box*(1-eps);
M_ref = M_max + 2;
opt.xi = pi*(M_max/L) / 15;
opt.P = 32;
opt.box = box;
opt.M = M_ref*box;
opt.oversampling = 1+sqrt(3);

% compute reference
opt.rc = inf;
ref = fse.direct_sum(x, f, box) - fse.real_sum(x, f, opt) - fse.self(f, opt);
ref_max = norm(ref(:), inf);
if do_3p
    ref_3p = fse.fourier_sum_3p(1:N, x, f, opt.xi, opt);
    ref_3p_max = norm(ref_3p(:), inf);
    ref_max = max(ref_max, ref_3p_max);
end
Mlist = 2:1:M_max;
[err_inf, err_rms, est, err_3p_rms] = deal(zeros(size(Mlist)));
tic
for i=1:numel(Mlist)
    this_opt = opt;
    this_opt.M = Mlist(i)*[1 1 1];
    
    fse_opt = setup_fse(this_opt);
    this_opt.R = fse_opt.R;
    
    uk = fse.fourier_sum(x, f, this_opt, fse.precomp(this_opt));    
    err = uk - ref;
    err = err;
    err_rms(i) = sqrt(1/N*sum(err(:).^2));
    est(i) = fse.fourier_estimate(f, this_opt);
    if do_3p
        est_3p(i) = fse.estimate_3p(f, this_opt);   
        u_3p = fourier_sum_3p(1:N, x, f, opt.xi, this_opt);
        err_3p = u_3p - ref_3p;
        err_3p_rms(i) = sqrt(1/N*sum(err_3p(:).^2));    
    end
end
toc

K = Mlist*pi/L;
xi = opt.xi;
Kxi = K / xi;

output.Kxi = Kxi;
output.err_rms_rel = err_rms / ref_max;
output.est_rel = est / ref_max;
output.xi = xi;

clf
semilogy(Kxi, err_rms / ref_max, '.', 'DisplayName','Measured')    
hold on
semilogy(Kxi, est / ref_max,'-', 'DisplayName','Estimate')
if do_3p
    semilogy(Kxi, est_3p / ref_max,'-', 'DisplayName','3P Estimate')
    semilogy(Kxi, err_3p_rms / ref_max, 'o', 'DisplayName','3P Measured')    
end
ylim([1e-16 1])
grid on
xlabel('$K/\xi$','interpreter','latex')
ylabel(sprintf('%s Fourier space RMS error (rel.)', kernel))
title(sprintf('N=%d, \\xi=%g, L=%d',N,opt.xi,L)) 
legend toggle
drawnow

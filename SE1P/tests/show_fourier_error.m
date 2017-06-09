function output = show_fourier_error(varargin)

p = inputParser;
p.addParameter('M',30);
p.addParameter('L',2);
p.addParameter('N',100);

p.parse(varargin{:});
input = p.Results;
M_max = input.M;
L = input.L;
N = input.N;

box = L*[1 1 1];
[x, f] = vector_system(N, box);

M_ref = M_max + 2;
opt.xi = pi*(M_max/L) / 12;
opt.box = box;
%opt.layers = 20;
opt.M = M_ref*box(3);
opt.P = 32;
opt.sl = 4;
opt.nl = input.M-2;
opt.s0 = 2.5;

% compute reference
ref = se1p_fourier_space(x, f, opt);
ref_max = rms(ref);

Mlist = 4:1:M_max;
[err_inf, err_rms, est] = deal(zeros(size(Mlist)));
tic
for i=1:numel(Mlist)
    this_opt = opt;
    this_opt.M = Mlist(i);
    uk = se1p_fourier_space(x, f, this_opt);
    err = uk - ref;
    err_rms(i) = rms(err);
    est(i) = fourier_estimate(f, this_opt);
end
toc
ref_max=1;
K = Mlist*pi/L;
xi = opt.xi;
Kxi = K / xi;

output.Kxi = Kxi;
output.err_rms_rel = err_rms / ref_max;
output.est_rel = est / ref_max;
output.xi = xi;

clf
semilogy(Kxi, err_rms / ref_max, 'b.', 'DisplayName','Measured')    
hold on
semilogy(Kxi, est / ref_max,'r-', 'DisplayName','Estimate')
ylim([1e-16 1])
grid on
xlabel('$K/\xi$','interpreter','latex')
ylabel('Fourier space RMS error (rel.)')
title(sprintf('N=%d, \\xi=%g, L=%d',N,opt.xi,L)) 
legend toggle
drawnow


function e = fourier_estimate(f, opt);
Q = sum(f(:).^2);
xi = opt.xi;
L = opt.box(1);
K = opt.M(1)/2;


e = xi*pi^(-2)*K^(-3/2)*sqrt(Q)*exp(-(pi*K/(xi*L))^2);
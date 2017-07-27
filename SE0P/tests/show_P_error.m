function output = show_P_error(kernel, varargin)

p = inputParser;
p.addParameter('Kxi',20);
p.addParameter('M',20);
p.addParameter('L',1);
p.addParameter('N',100);
p.addParameter('no_extra_support',false);
p.addParameter('oversampling', 1+sqrt(3));
p.addParameter('Plist',  [8 16 24]);
p.addParameter('M_min',2);

p.parse(varargin{:});
input = p.Results;
M_max = input.M;
Kxi_max = input.Kxi;
L = input.L;
N = input.N;

fse = fsewald(kernel);

opt.oversample_all = false;
opt.no_extra_support = input.no_extra_support;

Plist = input.Plist;
Mlist = input.M_min : 2 : M_max;

rng(1);
box = L*[1 1 1];
[x f] = fse.generator(N, box);
x(1,:) = 0; x(2,:) = box*(1-eps);
M_ref = M_max + 2;
opt.xi = pi*(M_max/L) / Kxi_max;
opt.P = 32;
opt.box = box;
opt.M = M_ref*box;
opt.oversampling = input.oversampling;

% compute reference
opt.rc = inf;
full_ref = fse.direct_sum(x, f, box);
full_ref_max = norm(full_ref(:), inf);
full_ref_rms = sqrt(1/N*sum(full_ref(:).^2));


ref = full_ref - fse.real_sum(x, f, opt) - fse.self(f, opt);
ref_max = norm(ref(:), inf);

% Iterate
[err_rms] = deal(zeros(numel(Plist),numel(Mlist)));
tic
for j=1:numel(Plist)
    opt.P = Plist(j);
    for i=1:numel(Mlist)
        this_opt = opt;
        this_opt.M = Mlist(i)*[1 1 1];
                
        uk = fse.fourier_sum(x, f, this_opt, fse.precomp(this_opt));    
        err = uk - ref;
        err = err;
        err_rms(j,i) = sqrt(1/N*sum(err(:).^2));
    end
end
toc

K = Mlist*pi/L;
xi = opt.xi;
Kxi = K / xi;

clf
styles = {'.-k', 'o-k', '+-k', '*-k', 's-k'};
for j=1:numel(Plist)
    semilogy(Kxi, err_rms(j,:) / ref_max, styles{mod(j-1, numel(styles))+1}, ...
             'DisplayName', sprintf('P=%d',Plist(j)))    
    hold on
end

if opt.oversample_all
    pstr = ', full oversampling';
else
    pstr = ', precomp oversampling';
end

ylim([1e-16 10])
grid on
xlabel('$\pi M / \xi L$','interpreter','latex')
ylabel(sprintf('%s Fourier space RMS error (rel.)', kernel))
title(sprintf('N=%d, \\xi=%g, L=%d, s_f=%g%s',N,opt.xi,L,opt.oversampling,pstr)) 
h = legend('toggle');
set(h, 'Location','SouthWest')
drawnow

output.Kxi = Kxi;
output.err_rms = err_rms;
output.ref_max = ref_max;
output.full_ref_max = full_ref_max;
output.full_ref_rms = full_ref_rms;
output.xi = xi;
output.Plist = Plist;
output.Mlist = Mlist;
function fse_opt = setup_fse(opt)
% Setup Free-Space Ewald grid sizes
% inner = domain containing all sources and targets
% extended = FGG grid (to cancel wrap effects)
% padded = 2x FFT grid (for aperiodic convolution)
% oversampled = FFT grid for truncated Green's function

inner_M = opt.M;
inner_box = opt.box;
h = opt.box(1) / opt.M(1);
popt = parse_params(opt);
P = popt.P;
m = popt.m;
xi = popt.xi;
eta = popt.eta;
% Old delta (should cover grid Gaussian)
delta_old = h*P/2;
% New delta (should cover remainder Gaussian)
if popt.eta < 1
    deltaNew = sqrt((1-eta)/2)*m/xi;
else
    deltaNew = delta_old;
end
% Take max
delta = max(deltaNew, delta_old);

if isfield(opt, 'no_extra_support') && opt.no_extra_support == true;
    delta = delta_old;
end
%fprintf('P=%d, eta=%.2f, delta1=%.2f, delta2=%.2f, delta=%.2f\n', ...
%        P, eta, delta_old, deltaNew, delta);

EVEN_GRIDS = true;
if EVEN_GRIDS
    deltaM = 2*ceil(delta / h);
else
    deltaM = ceil(2*delta / h);
end
delta = h*deltaM/2;

fse_opt.extended_box = inner_box + 2*delta;
fse_opt.extended_M = inner_M + deltaM;

fse_opt.padded_M = fse_opt.extended_M * 2;
fse_opt.padded_box = fse_opt.extended_box * 2;

% Ensure integer oversampling rate
if EVEN_GRIDS 
    overM = ceil(opt.oversampling * fse_opt.extended_M/2)*2;
else
    overM = ceil(opt.oversampling * fse_opt.extended_M);
end
actual_oversampling = overM / fse_opt.extended_M;
if abs(actual_oversampling - opt.oversampling) > 10*eps(actual_oversampling)
    warning('FSE:OversamplingIncreased',...
            'Oversampling rate increased to %g achieve integer grid', ...
            actual_oversampling);
end
fse_opt.oversampled_M = overM; 
fse_opt.oversampled_box = fse_opt.extended_box * actual_oversampling;

fse_opt.R = norm(fse_opt.extended_box);
fse_opt.delta = delta;

if isfield(opt, 'oversample_all') && opt.oversample_all == true;
    fse_opt.padded_M = fse_opt.oversampled_M;
    fse_opt.padded_box = fse_opt.oversampled_box;
end

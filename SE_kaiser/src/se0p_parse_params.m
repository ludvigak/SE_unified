function fse_opt = se0p_parse_params(opt)
% Setup Free-Space Ewald grid sizes
% inner = domain containing all sources and targets
% extended = FGG grid (to cancel wrap effects)
% padded = 2x FFT grid (for aperiodic convolution)
% oversampled = FFT grid for truncated Green's function

fse_opt = opt;
inner_M = opt.M;
inner_box = opt.box;
h = opt.box(1) / opt.M(1);

P = opt.P;
xi = opt.xi;

fse_opt.L = opt.box(1);
fse_opt.h = h;
fse_opt.p_half = (mod(fse_opt.P,2)==0)*fse_opt.P/2+(mod(fse_opt.P,2)~=0)*(fse_opt.P-1)/2;
fse_opt.p_half = fse_opt.P/2;

if(isfield(opt,'beta'))
    fse_opt.beta=opt.beta*P; 
else 
    fse_opt.beta=2.5*P; 
end

fse_opt.oversampling = opt.s;
delta = h*P/2;

EVEN_GRIDS = true;
if EVEN_GRIDS
    deltaM = 2*ceil(delta / h);
else
    deltaM = ceil(2*delta / h);
end

delta = h*deltaM/2;
% should cover the gaussian factor
delta = max(delta, sqrt(opt.beta*2)/xi);

fse_opt.extended_box = inner_box + 2*delta;
fse_opt.extended_M = inner_M + deltaM;

fse_opt.padded_M = fse_opt.extended_M * 2;
fse_opt.padded_box = fse_opt.extended_box * 2;

% Ensure integer oversampling rate
if EVEN_GRIDS 
    overM = ceil(fse_opt.oversampling * fse_opt.extended_M/2)*2;
else
    overM = ceil(fse_opt.oversampling * fse_opt.extended_M);
end
actual_oversampling = overM / fse_opt.extended_M;

fse_opt.oversampled_M = overM;
fse_opt.oversampled_box = fse_opt.extended_box * actual_oversampling;

fse_opt.R = norm(fse_opt.extended_box);
fse_opt.delta = delta;

if isfield(opt, 'oversample_all') && opt.oversample_all == true;
    fse_opt.padded_M = fse_opt.oversampled_M;
    fse_opt.padded_box = fse_opt.oversampled_box;
end
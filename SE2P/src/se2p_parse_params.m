function popt = se2p_parse_params(opt)
 
% check that we have all mandatory options
assert(isfield(opt,'M'))
assert(isfield(opt,'P'))
assert(isfield(opt,'box'))

% copy all params
popt = opt;

% step size
popt.L = opt.box(1);
popt.h = popt.L/popt.M;

% Gaussian
popt.m = .94*sqrt(pi*popt.P);
w = popt.h*popt.P/2;
eta = (2*w*popt.xi/popt.m)^2;
popt.c = 2*popt.xi^2/eta;

popt.w = popt.h*popt.P/2;
popt.Mz = ceil(popt.box(3)/popt.h)+popt.P;
popt.Lz = popt.box(3)+2*popt.w;

% Even grids if even grids are in the z direction
if(mod(popt.M,2)==0)
    popt.Mz = 2*ceil(popt.Mz/2);
end

% h should be the same in all directions
assert((popt.Lz/popt.Mz-popt.h)<eps)

% sampling factor (oversampling)
if( isfield(opt,'s')), popt.s = opt.s; else popt.s=1; end;

popt.R = popt.Lz;

% increase s and such that FFTN has integer size vectors.
popt.s = ceil(popt.s*popt.Mz)/popt.Mz;


% collect
popt.PH = popt.P/2;
wbox = [0 popt.L; 0 popt.L; -popt.w  -popt.w+popt.Lz];
popt.a = wbox(3,1);
popt.eta = (2*popt.w*popt.xi/popt.m)^2;
popt.c = 2*popt.xi^2/popt.eta;

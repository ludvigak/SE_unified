function opt = se3p_parse_params(opt)
assert(isfield(opt,'M'));
assert(isfield(opt,'xi'));
assert(isfield(opt,'P'))
assert(isfield(opt,'box'))

h = opt.box(1)/opt.M(1);
assert(abs( h-opt.box(2)/opt.M(2))<eps);
assert(abs( h-opt.box(3)/opt.M(3))<eps);

opt.h = h;
opt.w = h*(opt.P-1)/2;
opt.m = .98*sqrt(pi*opt.P);

%wbox = [-p.w  -popt.w+popt.Lx;-popt.w  -popt.w+popt.Ly; 0 popt.L];
opt.eta = (2*opt.w*opt.xi/opt.m)^2;
opt.c = 2*opt.xi^2/opt.eta;
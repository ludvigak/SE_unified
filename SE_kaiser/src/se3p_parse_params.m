function opt = se3p_parse_params(opt)
assert(isfield(opt,'M'));
assert(isfield(opt,'P'));
assert(isfield(opt,'box'));

opt.h = opt.box(1)/opt.M(1);
assert(abs( opt.h-opt.box(2)/opt.M(2))<eps)
assert(abs( opt.h-opt.box(3)/opt.M(3))<eps)
opt.m = 1.71*sqrt(opt.P);
opt.p_half = (mod(opt.P,2)==0)*opt.P/2+(mod(opt.P,2)~=0)*(opt.P-1)/2;
opt.w = opt.h*opt.P/2;

if(isfield(opt,'beta')), opt.beta=opt.beta*opt.p_half; else opt.beta=4.74*opt.p_half; end
opt.eta = (2*opt.w*opt.xi/opt.m)^2;
opt.c = 2*opt.xi^2/opt.eta;
function p = parse_params(opt)

% check that we have all mandatory options
assert(isfield(opt,'M'))
assert(isfield(opt,'P'))
assert(isfield(opt,'box'))

% verify all assumptions on parameters

% step size
L = opt.box(1);
h = L/opt.M(1);
assert(abs(opt.box(2)/opt.M(2) - opt.box(1)/opt.M(1)) < eps)
assert(abs(opt.box(3)/opt.M(3) - opt.box(2)/opt.M(2)) < eps)

% Gaussian
P = opt.P;
if( isfield(opt,'m')), m = opt.m; else m = 0.9*sqrt(pi*P); end;
w = h*(P-1)/2;

% collect
p.M=opt.M;
p.P = P;
p.w = w;
p.m = m;
p.box = opt.box;
p.L = L;
p.h = h;

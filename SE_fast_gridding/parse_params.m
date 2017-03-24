function p = parse_params(opt)
% Parse parameter struct for fast Gaussian gridding
%
% Mandatory parameters:
%
% :param opt.M: grid size (1 x 3)
% :param opt.P: Gaussian width
% :param opt.box: box size (1 x 3)

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
if( isfield(opt,'m'))
    m = opt.m; 
else 
    m = 0.9*sqrt(pi*P); 
end;
w = h*(P-1)/2;

if isfield(opt, 'xi')
    opt.eta = (2*w*opt.xi/m)^2;
    opt.c = 2*opt.xi^2/opt.eta;
end

% External eval points
if isfield(opt,'eval_x') && numel(opt.eval_x)
    eval_external = 1;
    eval_x = opt.eval_x;
else
    eval_external = 0;
    eval_x = [];
end


% Default is to keep params
p = opt;
p
% collect
p.P = P;
p.w = w;
p.m = m;
p.L = L;
p.h = h;
p.eval_external = eval_external;
p.eval_x = eval_x;
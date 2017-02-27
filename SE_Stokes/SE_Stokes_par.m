function u  = SE_Stokes_par(eval_idx,x,f,xi,opt)
% parfor:ed version

verb = false;

% parameters and constants
opt = parse_params(opt);
[w m M P] = unpack_params(opt);
eta = (2*w*xi/m)^2;
opt.c = 2*xi^2/eta;

[G{1:3}]=deal(complex(zeros(M)));
parfor i=1:3
    % to grid, transform and shift
    G{i} = fftshift( fftn( SE_fg_grid_mex(x,f(:,i), opt) ) );
    if isreal(G{i})
        cprintf(verb, 'Forcing G{%d} complex.\n',i);
        G{i}=complex(G{i});
    end
end
cprintf(verb, 'M = [%d %d %d] P = %d m=%d w=%f\n',M,P,m,w);
cprintf(verb, 'eta = %f\t a=%f\n', eta, pi^2/opt.c);

% multiply with modified greens function
[G{1:3}] = stokeslet_fast_k_scaling(G{1},G{2},G{3},xi,opt.box,eta);

parfor i=1:3
    u(:,i) = SE_fg_int_mex(x(eval_idx,:), real( ifftn( ifftshift( G{i} ))) ,opt);
end

% ------------------------------------------------------------------------------
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


% ------------------------------------------------------------------------------
function [w m M P] = unpack_params(opt)
w = opt.w;
m = opt.m;
M = opt.M;
P = opt.P;

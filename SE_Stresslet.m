function [u stats]  = SE_Stresslet(eval_idx,x,f,n,xi,opt,varargin)

addpath('bin')

verb = true;

% parameters and constants
opt = parse_params(opt);
box = opt.box;
[w m M P] = unpack_params(opt);
eta = (2*w*xi/m)^2;
opt.c = 2*xi^2/eta;

if nargin == 6
    static_fgg=false;
    sdat = [];
elseif nargin == 7
    static_fgg=true;
    sdat=varargin{1};
    x = x(sdat.perm,:);
    f = f(sdat.perm,:);
    n = n(sdat.perm,:);
end

% Grid
G=zeros([9 M]);

% Indices for outer product
n_idx = [1 2 3 1 2 3 1 2 3];
f_idx = [1 1 1 2 2 2 3 3 3];

tic;
parfor i=1:9
    % to grid, transform and shift
    i1 = n_idx(i);
    i2 = f_idx(i);
    S = n(:,i1).*f(:,i2);
    if static_fgg
        G( i, :, :, :) = fftshift( fftn( ...
                            SE_fg_grid_split_mex(x,S,opt, ...
                              sdat.zs,sdat.zx,sdat.zy,sdat.zz,sdat.idx ) ...
                         ) );
    else
        G( i, :, :, :) = fftshift( fftn( SE_fg_grid_mex(x,S,opt) ) );
    end
end
stats.wtime_grid = toc();

cprintf(verb, 'M = [%d %d %d] P = %d m=%d w=%f\n',M,P,m,w);
cprintf(verb, 'eta = %f\t a=%f\n', eta, pi^2/opt.c);

% Do scaling
tic;
[H{1:3}] = stresslet_fast_k_scaling(G,xi,box,eta); % Mex
% [H{1:3}] = stresslet_k_scaling(G,xi,box,eta); % Matlab
stats.wtime_scale = toc();

if static_fgg
    u = zeros(size(x));
    x = x(sdat.iperm,:);
else
    u = zeros(length(eval_idx),3);
end

% Back transform and gather
tic;
parfor i=1:3
    if static_fgg
        u(:,i) = SE_fg_int_split_mex(x,...
                                    real( ifftn( ifftshift( H{i} ))),opt,...
                                    sdat.zs,sdat.zx,sdat.zy,sdat.zz,sdat.idx);
    else
        u(:,i) = SE_fg_int_mex(x(eval_idx,:), ...
                                    real( ifftn( ifftshift( H{i} ))) ,opt);
    end
end
stats.wtime_int = toc();

if static_fgg
     u = u(sdat.iperm(eval_idx),:);
end

stats.wtime=[stats.wtime_grid stats.wtime_scale stats.wtime_int];
stats.wtime_total = sum(stats.wtime);


% ------------------------------------------------------------------------------
function [w m M P] = unpack_params(opt)
w = opt.w;
m = opt.m;
M = opt.M;
P = opt.P;
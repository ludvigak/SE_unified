function [u stats]  = SE_Stresslet(eval_idx,x,f,n,xi,opt,varargin)

ttic = tic;

verb = false;

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
G=complex(zeros([9 M]));

% Indices for outer product
n_idx = [1 2 3 1 2 3 1 2 3];
f_idx = [1 1 1 2 2 2 3 3 3];

% Save timings for gridding+FFT separately
gftic = tic;
[gtime ftime] = deal(zeros(9,1));

% Parfor speeds up split gridding code a little bit,
% but actually slows down regular gridding
parfor i=1:9
    % to grid, transform and shift
    i1 = n_idx(i);
    i2 = f_idx(i);
    S = n(:,i1).*f(:,i2);
    gtic = tic;
    if static_fgg
        F = SE_fg_grid_split_mex(x,S,opt, ...
                                sdat.zs,sdat.zx,sdat.zy,sdat.zz,sdat.idx );
    else
        F = SE_fg_grid_mex(x,S,opt);
    end
    gtime(i) = toc(gtic);
    ftic = tic;
    G( i, :, :, :) = fftshift( fftn( F ) );
    ftime(i) = toc(ftic);
end

wgridfft = toc(gftic); % Total time spent in loop
gtime = sum(gtime); % Total time spent on gridding (over all threads)
ftime = sum(ftime); % Total time spent on FFT (over all threads)
% Assume that all thread time was spent on gridding + FFT
stats.wtime_grid = wgridfft*gtime/(gtime+ftime);
stats.wtime_fft = wgridfft*ftime/(gtime+ftime);


cprintf(verb, 'M = [%d %d %d] P = %d m=%d w=%f\n',M,P,m,w);
cprintf(verb, 'eta = %f\t a=%f\n', eta, pi^2/opt.c);

if isreal(G)
    cprintf(verb,'Forcing G complex.\n');
    G = complex(G);
end

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
iftic = tic;
[itime ftime] = deal(zeros(3,1));

% Integration code parallelized using OpenMP in MEX
for i=1:3
    ftic = tic;
    F = real( ifftn( ifftshift( H{i} )));
    ftime(i) = toc(ftic);
    itic = tic;
    if static_fgg
        u(:,i) = SE_fg_int_split_mex(x,F,opt,...
                                    sdat.zs,sdat.zx,sdat.zy,sdat.zz,sdat.idx);
    else
        u(:,i) = SE_fg_int_mex(x(eval_idx,:), F ,opt);
    end
    itime(i) = toc(itic);
end

wintfft = toc(iftic);
itime = sum(itime);
ftime = sum(ftime);
stats.wtime_int = wintfft*itime/(itime+ftime);
stats.wtime_fft = wintfft*ftime/(itime+ftime) + stats.wtime_fft;

if static_fgg
     u = u(sdat.iperm(eval_idx),:);
end

stats.wtime_total = toc(ttic);


% ------------------------------------------------------------------------------
function [w m M P] = unpack_params(opt)
w = opt.w;
m = opt.m;
M = opt.M;
P = opt.P;
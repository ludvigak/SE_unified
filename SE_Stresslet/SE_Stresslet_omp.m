function [u stats]  = SE_Stresslet_omp(eval_idx,x,f,n,xi,opt,varargin)

ttic = tic;

verb = false;

% parameters and constants
opt = parse_params(opt);
box = opt.box;
[w m M P] = unpack_params(opt);
eta = (2*w*xi/m)^2;
opt.c = 2*xi^2/eta;

x = recenter_points(x, box);

if nargin == 6
    static_fgg=false;
    sdat = [];
elseif nargin == 7
    %warning('TBD:SE_STATIC','SE_static not passed to gridding');
    static_fgg=true;
    sdat=varargin{1};
    x = x(sdat.perm,:);
    f = f(sdat.perm,:);
    n = n(sdat.perm,:);
end

% Grid
gtic = tic;
G = stresslet_fg_grid_mex(x,n,f,opt);
stats.wtime_grid = toc(gtic);

% FFT
ftic = tic;
for i=1:9
    G(:,:,:,i) = fftshift( fftn( G(:,:,:,i) ) );
end
stats.wtime_fft = toc(ftic); % Total time spent on FFT

% Shuffle
shtic = tic();
G = permute(G, [4 1 2 3]);
stats.wtime_shuffle = toc(shtic);


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

if opt.eval_external
    u = zeros(size(opt.eval_x));
    opt.eval_x = recenter_points(opt.eval_x, box);
elseif static_fgg
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
    if opt.eval_external
        u(:,i) = SE_fg_int_mex(opt.eval_x, F ,opt);
    elseif static_fgg
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

if static_fgg && opt.eval_external==0
     u = u(sdat.iperm(eval_idx),:);
end

stats.wtime_total = toc(ttic);


% ------------------------------------------------------------------------------
function [w m M P] = unpack_params(opt)
w = opt.w;
m = opt.m;
M = opt.M;
P = opt.P;
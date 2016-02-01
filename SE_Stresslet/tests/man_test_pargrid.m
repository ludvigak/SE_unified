clear

% test effect of parallelizing gridding with parfor (non-SSE version)

box = [1 1 1];

% grid
opt.M = 100*box;
opt.box = box;
opt.L = opt.box(1);
opt.h = opt.L/opt.M(1);

xi = 4;        % ewald parameter
P = 16;
opt.P = P;
opt.m = 0.9*sqrt(pi*P);
opt.w = opt.h*(P-1)/2;
eta = (2*opt.w*xi/opt.m)^2;
opt.c = 2*xi^2/eta;

% Grid
G=zeros([9 opt.M]);

% Indices for outer product
n_idx = [1 2 3 1 2 3 1 2 3];
f_idx = [1 1 1 2 2 2 3 3 3];

Nlist = 10*2.^[4 5 6 7 8 9 10 11 12];
Tlist = [];
tests = 2;

for N=Nlist
    % system
    [x f n] = generate_state(N,box);	

    this_partime = [];
    for t=1:tests
        tic;
        parfor i=1:9
            % to grid, transform and shift
            i1 = n_idx(i);
            i2 = f_idx(i);
            S = n(:,i1).*f(:,i2);
            G( i, :, :, :) = fftshift( fftn( SE_fg_grid_mex(x,S,opt) ) );
        end
        this_partime(end+1) = toc;
    end
    
    this_sertime = []; 
    for t=1:tests  
        tic;
        for i=1:9
            % to grid, transform and shift
            i1 = n_idx(i);
            i2 = f_idx(i);
            S = n(:,i1).*f(:,i2);
            G( i, :, :, :) = fftshift( fftn( SE_fg_grid_mex(x,S,opt) ) );
        end
        this_sertime(end+1) = toc;
    end
    
    Tlist = [Tlist; mean(this_sertime) mean(this_partime)];
end

speedup = Tlist(:,1)./Tlist(:,2);
disp('-------------------------')
disp('    N    ser   par   spup')  
fprintf('% 7d  %.2f  %.2f  %.2f\n',[Nlist' Tlist speedup]')

fprintf('\nRun on %d cores.\n',matlabpool('size'));


%% M = 31
%**** Results for various compiles

% -msse2
%-------------------------
%    N    ser   par   spup
%    160  0.07  0.07  1.08
%    320  0.08  0.07  1.16
%    640  0.09  0.07  1.25
%   1280  0.12  0.08  1.51
%   2560  0.19  0.11  1.75
%   5120  0.33  0.14  2.34
%  10240  0.57  0.24  2.44
%  20480  1.09  0.42  2.60
%  40960  2.11  0.79  2.66
%  81920  4.16  1.54  2.70
%
% Run on 4 cores.
 
% -xSSE4.1 (compiled on Ferlin)
% -------------------------
%     N    ser   par   spup
%     160  0.08  0.08  1.00
%     320  0.08  0.07  1.11
%     640  0.09  0.08  1.24
%    1280  0.13  0.09  1.39
%    2560  0.19  0.11  1.70
%    5120  0.31  0.14  2.18
%   10240  0.56  0.24  2.32
%   20480  1.05  0.41  2.54
%   40960  2.04  0.78  2.61
%   81920  4.01  1.51  2.66
% 
% Run on 4 cores.

%% M = 100
% -------------------------
%     N    ser   par   spup
%     160  1.24  1.47  0.85
%     320  1.41  1.50  0.94
%     640  1.48  1.52  0.97
%    1280  1.62  1.66  0.98
%    2560  1.87  1.84  1.01
%    5120  2.36  2.31  1.02
%   10240  3.26  3.18  1.02
%   20480  5.18  5.06  1.02
%   40960  9.14  8.59  1.06
% 
% Run on 4 cores.

% No gain on large grids!
% Probably too many cache misses when points are not sorted



function show_P_error
    global NOPRE N

restoredefaultpath
init

rng(1)
box = [1 1 1];   % domain
N = 10;          % number of charged particles
M0 = 28;

opt.M = M0*box;
opt.P = 20;
opt.xi = pi*M0 / 12;
opt.rc = 6 / opt.xi;
opt.box = box;
opt.layers = (opt.M(1)-1)/2;
opt.beta = 1.5*pi*.98^2;

% repetition
rep = 3;

% exclude precomputation
NOPRE = true;

% charge-neutral system
[x q] = vector_system(N,box);

Plist = 4:2:18;
[kaiser_3p gaussian_3p kaiser_2p gaussian_2p ...
 kaiser_1p gaussian_1p kaiser_0p gaussian_0p] = deal([]);
[tg_3p tk_2p tg_2p tk_1p tg_1p tk_0p tg_0p] = deal([]);

kaiser = zeros(4,numel(Plist));
tk     = zeros(4,numel(Plist));
gaussian= zeros(4,numel(Plist));
tg     = zeros(4,numel(Plist));

%% --------------------------- Iterate Kaiser ------------------------------
%% 3d-periodic
opt_3p = opt;
ref3p = compute_direct_3p(x,q,opt_3p);
for k=1:numel(Plist)
    opt_3p.P = Plist(k);
    [~, tk(4,k) kaiser(4,k)] = call_kaiser_3p(x, q, opt_3p, ref3p, rep);
end

%% 2d-periodic
opt_2p = set2p_params(opt);
ref2p = compute_direct_2p(x,q,opt_2p);
for k=1:numel(Plist)    
    opt_2p.P = Plist(k);
    [~, tk(3,k) kaiser(3,k)] = call_kaiser_2p(x, q, opt_2p, ref2p, rep);
end

%% 1d-periodic
opt_1p = set1p_params(opt);
ref1p = compute_direct_1p(x,q,opt_1p);
for k=1:numel(Plist)
    opt_1p.P = Plist(k);
    [~, tk(2,k) kaiser(2,k)] = call_kaiser_1p(x, q, opt_1p, ref1p, rep);
end

%% 0d-periodic
opt_0p = set0p_params(opt);
ref0p = compute_direct_0p(x,q,opt_0p);
for k=1:numel(Plist)    
    opt_0p.P = Plist(k);
    [~, tk(1,k) kaiser(1,k)] = call_kaiser_0p(x, q, opt_0p, ref0p, rep);
end


%% ------------------------------ Iterate gaussian --------------------------
%% 3d-periodic
opt_3p = opt;
addpath('../SE')
for k=1:numel(Plist)
    opt_3p.P = Plist(k);
    [~, tg(4,k) gaussian(4,k)] = call_gaussian_3p(x, q, opt_3p, ref3p, rep);
end

%% 2d-periodic
opt_2p = set2p_params(opt);
addpath('../SE2P/src')
for k=1:numel(Plist)    
    opt_2p.P = Plist(k);
    [~, tg(3,k) gaussian(3,k)] = call_gaussian_2p(x, q, opt_2p, ref2p, rep);
end

%% 1d-periodic
opt_1p = set1p_params(opt);
addpath('../SE1P/src')
for k=1:numel(Plist)
    opt_1p.P = Plist(k);
    [~, tg(2,k) gaussian(2,k)] = call_gaussian_1p(x, q, opt_1p, ref1p, rep);
end

%% 0d-periodic
opt_0p = set0p_params(opt);
addpath('../SE0P/src')
for k=1:numel(Plist)    
    opt_0p.P = Plist(k);
    [~, tg(1,k) gaussian(1,k)] = call_gaussian_0p(x, q, opt_0p, ref0p, rep);
end

%% -------------------------------- PLOT ------------------------------

colk={'bo-','bs-','bd-','b>-'};
colg={'ro-','rs-','rd-','r>-'};
      
sfigure(1);
hold on
for k=1:4
    plot(Plist,kaiser(k,:),colk{k},'DisplayName',sprintf('Kiaser %dp',k-1))
    plot(Plist,gaussian(k,:),colg{k},'DisplayName',sprintf('Gaussian %dp',k-1),...
	'markerFaceColor','r')
end
set(gca,'yscale','log')
xlabel('P')
ylabel('e_{rms} (rel.)')
h = legend('toggle');
set(h, 'Location','SouthWest')

sfigure(2);
hold on
for k=1:4
    plot(Plist,tk(k,:),colk{k},'DisplayName',sprintf('Kiaser %dp',k-1))
    plot(Plist,tg(k,:),colg{k},'DisplayName',sprintf('Gaussian %dp',k-1),...
	'markerFaceColor','r')
end
set(gca,'yscale','log')
xlabel('P')
ylabel('time [s]')
h = legend('toggle');
set(h, 'Location','SouthWest')

sfigure(3);
hold on
for k=1:4
    plot(kaiser(k,:),tk(k,:),colk{k},'DisplayName',sprintf('Kiaser %dp',k-1))
    plot(gaussian(k,:),tg(k,:),colg{k},'DisplayName',sprintf('Gaussian %dp',k-1),...
	'markerFaceColor','r')
end
set(gca,'xscale','log')
xlabel('time [s]')
xlabel('e_{rms} (rel.)')
h = legend('toggle');
set(h, 'Location','SouthWest')


end


%%% ---------------------------- SET params ------------------------------
function popt = set2p_params(opt)
    popt = opt;
    popt.M = opt.M(1)-2;
    popt.P = 20;
    popt.s = 3.5;
    popt.s0 = 2;
    popt.n = 6;    
end

function popt = set1p_params(opt)
    popt = opt;
    popt.M = opt.M(1);
    popt.P = 20;
    popt.s = 3.5;
    popt.s0 = 2.5;
    popt.n = 8;    
end

function popt = set0p_params(opt)
    popt = opt;
    popt.M = opt.M;
    popt.s = 2.7;
end

%%% ---------------------------- COMPUTE DIRECT ------------------------------
%% 3d-periodic direct
function ref3p = compute_direct_3p(x,q,opt)
    global N
% idx = 1:10;
%ref3p = SE3P_direct_fd_mex(idx,x,q,opt);
    ref3p = se3p_fourier_space_kaiser(1:N,x,q,opt);
end

%% 2d-periodic direct
function ref2p = compute_direct_2p(x,q,opt)
    global N
% idx = 1:10;
% ref2p = SE2P_direct_fd_mex(idx,x,q,opt) + SE2P_direct_k0_mex(idx,x,q,opt);
    ref2p = se2p_fourier_space_kaiser(1:N,x,q,opt);
end

%% 1d-periodic direct
function ref1p = compute_direct_1p(x,q,opt)
    global N
    %idx = 1:10;
 %ref1p = SE1P_direct_fd_mex(idx,x,q,opt) + SE1P_direct_k0_mex(idx,x,q,opt);
 ref1p = se1p_fourier_space_kaiser(1:N,x,q,opt);
end

%% 0d-periodic direct
function ref0p = compute_direct_0p(x,q,opt)
    ref0p = se0p_fourier_space_kaiser(x,q,opt);
end

%%% ---------------------------- CALL KAISER ------------------------------
function [u tk_3p kaiser_3p] = call_kaiser_3p(x, q, opt, ref3p, rep)
    global NOPRE N

    [u time] = se3p_fourier_space_kaiser(1:N,x,q,opt);
    kaiser_3p = rms(u-ref3p)/rms(ref3p);
    for r = 1:rep
        [u time] = se3p_fourier_space_kaiser(1:N,x,q,opt);
    end
    tk_3p = time.total-NOPRE*(time.pre+time.prefft);
end

function [u tk_2p kaiser_2p] = call_kaiser_2p(x, q, opt, ref2p, rep)
    global NOPRE N
    
    if(opt.P<8)
        d = 8;
    elseif(opt.P<12)
        d = 6;
    elseif(opt.P<16)
        d = 2;
    else
        d = 2;
    end
    opt.M = opt.M-d;
    [u time]= se2p_fourier_space_kaiser(1:N,x,q,opt);
    kaiser_2p   = rms(u-ref2p)/rms(ref2p);
    for r = 1:rep
        [u time]= se2p_fourier_space_kaiser(1:N,x,q,opt);
    end
    tk_2p = time.total-NOPRE*(time.pre+time.prefft);
end

function [u tk_1p kaiser_1p] = call_kaiser_1p(x, q, opt, ref1p, rep)
    global NOPRE N
    
    if(opt.P<14)
        d = 4;
    else
        d = 0;
    end
    opt.M = opt.M-d;
    [u time] = se1p_fourier_space_kaiser(1:N,x,q,opt);
    kaiser_1p= rms(u-ref1p)/rms(ref1p);
    for r = 1:rep
        [u time] = se1p_fourier_space_kaiser(1:N,x,q,opt);
    end
    tk_1p = time.total-NOPRE*(time.pre+time.prefft);
end

function [u tk_0p kaiser_0p] = call_kaiser_0p(x, q, opt, ref0p, rep)
    global NOPRE N

    if(opt.P<=10)
        d = 6;
    elseif(opt.P==12)
        d = 4;
    else
        d=2;
    end
    opt.M = opt.M-d;
    [u time] = se0p_fourier_space_kaiser(x,q,opt);
    kaiser_0p= rms(u-ref0p)/rms(ref0p);
    for r = 1:rep
        [u time] = se0p_fourier_space_kaiser(x,q,opt);
    end
    tk_0p = time.total-NOPRE*(time.pre+time.prefft);    
end

%%% ---------------------------- CALL GAUSSIAN ------------------------------
function [u tg_3p gaussian_3p] = ...
        call_gaussian_3p(x, q, opt, ref3p, rep)
    
    global NOPRE
    
    %    u    = spectral_ewald(1:N,x,q,opt.xi,opt);
    u   = se3p_fourier_space(x,q,opt);
    gaussian_3p   = rms(u-ref3p)/rms(ref3p);
    for r = 1:rep
        [u  time]  = se3p_fourier_space(x,q,opt);
    end
    tg_3p = time.total-NOPRE*(time.pre);
end
    
function [u tg_2p gaussian_2p] = call_gaussian_2p(x, q, opt, ref2p, rep)
    global NOPRE
    
    if(opt.P<16)
        d = 4;
    else
        d = 0;
    end
    opt.M = opt.M-d;
    [u time] = se2p_fourier_space(x,q,opt);
    gaussian_2p   = rms(u-ref2p)/rms(ref2p);
    for r = 1:rep
        [u time] = se2p_fourier_space(x,q,opt);
    end
    tg_2p = time.total-NOPRE*(time.pre);
end

function [u tg_1p gaussian_1p] = call_gaussian_1p(x, q, opt, ref1p, rep)
    global NOPRE
    opt.sl = opt.s;
    opt.nl=opt.n;
    [u time]    = se1p_fourier_space(x,q,opt);
    gaussian_1p  = rms(u-ref1p)/rms(ref1p);
    for r = 1:rep
        [u time]    = se1p_fourier_space(x,q,opt);
    end
    tg_1p = time.total-NOPRE*(time.pre);
end

function [u tg_0p gaussian_0p] = call_gaussian_0p(x, q, opt, ref0p, rep)
    global NOPRE
    fse_warnings('off')
    opt.oversampling = opt.s;
    pre = laplace_precomp(opt);
    u = laplace_fourier_space(x,q,opt,pre);
    gaussian_0p   = rms(u-ref0p)/rms(ref0p);
    for r = 1:rep
        t = tic;
        u  = laplace_fourier_space(x,q,opt,pre);
        time = toc(t);
    end
    tg_0p = time;
end
function demo_full
    global N

restoredefaultpath
init

rng(1)
box = [1 1 1];   % domain
N = 2e4;          % number of charged particles
M0 = 24;

opt.M = M0*box;
opt.P = 20;
opt.xi = pi*M0 / 12;
opt.rc = 6 / opt.xi;
opt.box = box;
opt.layers = (opt.M(1)-1)/2;
opt.beta = 2.5;

% charge-neutral system
[x q] = vector_system(N,box);

Pg = 24;
Pk = 16;

%% --------------------------- KAISER ------------------------------
%% 3d-periodic
opt_3p = opt;
disp('direct 3p')
ref3p = compute_direct_3p(x,q,opt_3p);
disp('Kaiser 3p')
opt_3p.P = Pk;
[~, kaiser(4)] = call_kaiser_3p(x, q, opt_3p, ref3p);


%% 2d-periodic
opt_2p = set2p_params(opt);
disp('direct 2p')
ref2p = compute_direct_2p(x,q,opt_2p);
disp('Kaiser 2p')
opt_2p.P = Pk;
[~, kaiser(3)] = call_kaiser_2p(x, q, opt_2p, ref2p);

%% 1d-periodic
opt_1p = set1p_params(opt);
disp('direct 1p')
ref1p = compute_direct_1p(x,q,opt_1p);
disp('Kaiser 1p')
opt_1p.P = Pk;
[~, kaiser(2)] = call_kaiser_1p(x, q, opt_1p, ref1p);

%% 0d-periodic
opt_0p = set0p_params(opt);
disp('direct 0p')
ref0p = compute_direct_0p(x,q,opt_0p);
disp('Kaiser 0p')
opt_0p.P = Pk;
[~, kaiser(1)] = call_kaiser_0p(x, q, opt_0p, ref0p);

%% ------------------------------ GAUSSIAN --------------------------
%% 3d-periodic
disp('gaussian 3p')
opt_3p = opt;
addpath('../SE')
opt_3p.P = Pg;
[~, gaussian(4)] = call_gaussian_3p(x, q, opt_3p, ref3p);

disp('gaussian 2p')
%% 2d-periodic
opt_2p = set2p_params(opt);
addpath('../SE2P/src')
opt_2p.P = Pg;
[~, gaussian(3)] = call_gaussian_2p(x, q, opt_2p, ref2p);

disp('gaussian 1p')
%% 1d-periodic
opt_1p = set1p_params(opt);
addpath('../SE1P/src')
opt_1p.P = Pg;
[~, gaussian(2)] = call_gaussian_1p(x, q, opt_1p, ref1p);

disp('gaussian 0p')
%% 0d-periodic
opt_0p = set0p_params(opt);
addpath('../SE0P/src')
opt_0p.P = Pg;
[~, gaussian(1)] = call_gaussian_0p(x, q, opt_0p, ref0p);

var={'Function','Three_P','Two_P','One_P', 'FreeSp'};
table({'KAISER','GAUSSIAN'}',[kaiser(4) gaussian(4)]',[kaiser(3) gaussian(3)]',...
      [kaiser(2) gaussian(2)]',[kaiser(1) gaussian(1)]','variablenames',var)

end

%%% ---------------------------- SET params ------------------------------
function popt = set2p_params(opt)
    popt = opt;
    popt.M = opt.M(1);
    popt.P = 20;
    popt.s = 4;
    popt.s0 = 2.5;
    popt.n = 8;
end

function popt = set1p_params(opt)
    popt = opt;
    popt.M = opt.M(1);
    popt.P = 20;
    popt.s = 4;
    popt.s0 = 2.7;
    popt.n = 8;
end

function popt = set0p_params(opt)
    popt = opt;
    popt.P = 20;
    popt.s = 2.8;
    popt.M = opt.M + 2;
end

%%% ---------------------------- COMPUTE DIRECT ------------------------------
%% 3d-periodic direct
function ref3p = compute_direct_3p(x,q,opt)
    global N
    ref3p = se3p_fourier_space_kaiser(1:N,x,q,opt);
end

%% 2d-periodic direct
function ref2p = compute_direct_2p(x,q,opt)
    global N
    ref2p = se2p_fourier_space_kaiser(1:N,x,q,opt);
end

%% 1d-periodic direct
function ref1p = compute_direct_1p(x,q,opt)
    global N
 ref1p = se1p_fourier_space_kaiser(1:N,x,q,opt);
end

%% 0d-periodic direct
function ref0p = compute_direct_0p(x,q,opt)
    ref0p = se0p_fourier_space_kaiser(x,q,opt);
end

%%% ---------------------------- CALL KAISER ------------------------------
function [u kaiser_3p] = call_kaiser_3p(x, q, opt, ref3p)
    global N
    u = se3p_fourier_space_kaiser(1:N,x,q,opt);
    kaiser_3p = rms(u-ref3p)/rms(ref3p);
end

function [u kaiser_2p] = call_kaiser_2p(x, q, opt, ref2p)
    global N
    
    u = se2p_fourier_space_kaiser(1:N,x,q,opt);
    kaiser_2p   = rms(u-ref2p)/rms(ref2p);
end

function [u kaiser_1p] = call_kaiser_1p(x, q, opt, ref1p)
    global N
    
    u = se1p_fourier_space_kaiser(1:N,x,q,opt);
    kaiser_1p= rms(u-ref1p)/rms(ref1p);
end

function [u kaiser_0p] = call_kaiser_0p(x, q, opt, ref0p)
    global N

    u = se0p_fourier_space_kaiser(x,q,opt);
    kaiser_0p= rms(u-ref0p)/rms(ref0p);
end

%%% ---------------------------- CALL GAUSSIAN ------------------------------
function [u gaussian_3p] = call_gaussian_3p(x, q, opt, ref3p)   
    u   = se3p_fourier_space(x,q,opt);
    gaussian_3p   = rms(u-ref3p)/rms(ref3p);
end
    
function [u gaussian_2p] = call_gaussian_2p(x, q, opt, ref2p)
    u = se2p_fourier_space(x,q,opt);
    gaussian_2p   = rms(u-ref2p)/rms(ref2p);
end

function [u gaussian_1p] = call_gaussian_1p(x, q, opt, ref1p)
    opt.sl = opt.s;
    opt.nl=opt.n;
    u     = se1p_fourier_space(x,q,opt);
    gaussian_1p  = rms(u-ref1p)/rms(ref1p);
end

function [u gaussian_0p] = call_gaussian_0p(x, q, opt, ref0p)
    fse_warnings('off')
    opt.oversampling = opt.s;
    pre = laplace_precomp(opt);
    u = laplace_fourier_space(x,q,opt,pre);
    gaussian_0p   = rms(u-ref0p)/rms(ref0p);
end

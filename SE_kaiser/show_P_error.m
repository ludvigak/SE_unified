clc; clear all

restoredefaultpath
init

rng(2)
box = [1 1 1];   % domain
N = 100;          % number of charged particles

M0 = 28;

opt.M = M0*box;
opt.xi = pi*M0 / 12;
opt.rc = 6 / opt.xi;
opt.box = box;

opt.layers = (opt.M(1)-1)/2;

opt.beta = 1.4*pi*.98^2;

% charge-neutral system
[x q] = vector_system(N,box);

idx = 1:min(N,50);
% compute FD Ewald sum
ref1p = SE1P_direct_fd_mex(idx,x,q,opt) + SE1P_direct_k0_mex(idx,x,q,opt);
ref2p = SE2P_direct_fd_mex(idx,x,q,opt) + SE2P_direct_k0_mex(idx,x,q,opt);
ref3p = SE3P_direct_fd_mex(idx,x,q,opt);

Plist = 4:2:18;
[kaiser_3p gaussian_3p kaiser_2p gaussian_2p kaiser_1p gaussian_1p] ...
    = deal([]);

for P=Plist
    opt.P = P;
    opt.M = M0*box;
    [u time] = se3p_fourier_space_kaiser(1:N,x,q,opt);
    kaiser_3p(end+1) = rms(u(idx)-ref3p)/rms(ref3p);

    if(P<8)
        d = 10;
    elseif(P<12)
        d = 8;
    elseif(P<16)
        d = 6;
    elseif(P==16)
        d = 4;
    else
        d = 2;
    end
    opt.M = M0*box(1)-d;
    opt.s = 3.5;
    opt.s0= 2.5;
    opt.n = 5;
    [u time]= se2p_fourier_space_kaiser(1:N,x,q,opt);
    kaiser_2p(end+1)   = rms(u(idx)-ref2p)/rms(ref2p);

    if(P<14)
        d = 4;
    else
        d = 2;
    end
    opt.M = M0*box(1)-d;
    opt.s = 3.5;
    opt.s0= 2.5;
    opt.n = 5;
    [u time] = se1p_fourier_space_kaiser(1:N,x,q,opt);
    kaiser_1p(end+1)= rms(u(idx)-ref1p)/rms(ref1p);
end

addpath('../SE')
for P=Plist
    opt.P = P;
    opt.M = M0*box;
    u    = spectral_ewald(1:N,x,q,opt.xi,opt);
    gaussian_3p(end+1)   = rms(u(idx)-ref3p)/rms(ref3p);
end

addpath('../SE2P/src')
for P=Plist
    opt.P = P;
    if(P<16)
        d = 4;
    else
        d = 0;
    end
    opt.M = M0*box(1)-d;
    u    = se2p_fourier_space(x,q,opt);
    gaussian_2p(end+1)   = rms(u(idx)-ref2p)/rms(ref2p);
end

addpath('../SE1P/src')
for P=Plist
    opt.P = P;
    opt.M = M0*box(1);
    opt.sl=opt.s;
    opt.nl=5;
    u      = se1p_fourier_space(x,q,opt);
    gaussian_1p(end+1)   = rms(u(idx)-ref1p)/rms(ref1p);
end

semilogy(Plist,kaiser_3p,'b.-',Plist,gaussian_3p,'ro-')
hold on
semilogy(Plist,kaiser_2p,'bs-',Plist,gaussian_2p,'rs-')
semilogy(Plist,kaiser_1p,'b^-',Plist,gaussian_1p,'r^-')
legend('Kaiser      3P','Gaussian 3P','Kaiser      2P', ...
       'Gaussian 2P',['Kaiser      1P'],'Gaussian 1P')
xlabel('P')
ylabel('e_{rms} (rel.)')

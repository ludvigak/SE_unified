% rms error in computation of the force as a function of P
% for different sl values and large nl.

clear
rng(1);

N = 200;
L = 1;
box = [L L L];
[x, f] = vector_system(N, box);

opt.M = 38;
opt.box = box;
opt.layers = 20;
opt.xi = 8;
opt.nl = 5;
opt.s0 = 2.5;

% direct
uf = SE1P_direct_fd_mex(1:N,x,f,opt);
u0 = SE1P_direct_k0_mex(1:N,x,f,opt);
ref = uf + u0;

sty = {'go-','b.-','rs-','ks-'};
Plist = 16:2:28;
figure(1)
clf
for sl = [1 2 3 4]
    opt.sl = sl;
    error = [];
    for p=Plist
        opt.P = p;    
        u   = se1p_fourier_space(x, f, opt);
        error(end+1) = rms(u-ref) / rms(ref);
    end
    name = sprintf('s_l=%d',sl);
    semilogy(Plist,error,sty{sl},'DisplayName',name)
    xlabel('P')
    ylabel('e_{rms} (rel.)')
    hold on
end
legend toggle

figure(2)
clf
opt.sl = 4;
for nl = [1 3 5]
    error = [];
    opt.nl = nl;
    for p=Plist
        opt.P = p;    
        u   = se1p_fourier_space(x, f, opt);
        error(end+1) = rms(u-ref) / rms(ref);
    end
    name = sprintf('n_l=%d',nl);
    semilogy(Plist,error,sty{(nl+1)/2+1},'DisplayName',name)
    xlabel('P')
    ylabel('e_{rms} (rel.)')
    hold on
end
legend toggle
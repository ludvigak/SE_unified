clear
rng(1);

N = 100;
L = 1;
box = [L L L];
[x, f] = vector_system(N, box);

M0 = 28; % Set M0 to an even number, the rest is auto

opt.M = M0*box(1);
opt.box = box;
opt.layers = 20;
opt.xi = pi*M0 / 12;
opt.nl = 4;
opt.sl = 3;
opt.s0 = 2.5;

% direct
uf = SE1P_direct_fd_mex(1:N,x,f,opt);
u0   = SE1P_direct_k0_mex(1:N,x,f,opt);
ref = uf + u0;

error = [];
Plist = 4:4:24;
for p=Plist
    opt.P = p;    
    u   = se1p_fourier_space(x, f, opt);
    error(end+1) = rms(u-ref) / rms(ref);
end
semilogy(Plist,error,'.-')
xlabel('P')
ylabel('e_{rms} (rel.)')
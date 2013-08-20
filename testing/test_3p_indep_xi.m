clear all, close all

N = 10;
opt.layers = 30;
opt.box = [1.1 1.32 1.192];

x = [opt.box(1)*rand(N,1) opt.box(2)*rand(N,1) opt.box(3)*rand(N,1)];
f = rand(N,3); f = f-repmat(mean(f),N,1);
idx = 1;

xi = linspace(1,10,5);

for j=1:length(xi)
    opt.xi = xi(j);
    u1 = SE3P_Stokes_direct_real_mex(idx,x,f,opt);
    u2 = SE3P_Stokes_direct_fd_mex(idx,x,f,opt);
    u3 = SE3P_Stokes_direct_self_mex(idx,x,f,opt);
    u(j,:) = u1+u2+u3;
end
u
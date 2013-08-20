clear all, close all

N = 10;
opt.layers = 30;
opt.box = [1.1 1.32 1.192];

x = [opt.box(1)*rand(N,1) opt.box(2)*rand(N,1) opt.box(3)*rand(N,1)];
q = rand(N,1); q = q-mean(q);
idx = 1;

xi = linspace(1,10,5);

for j=1:length(xi)
    opt.xi = xi(j);
    phi1 = SE3P_direct_real_mex(idx,x,q,opt);
    phi2 = SE3P_direct_fd_mex(idx,x,q,opt);
    phi3 = SE3P_direct_self_mex(idx,x,q,opt);
    phi(j) = phi1+phi2+phi3;
end
phi(:)
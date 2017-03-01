clear all, close all

N = 10;
opt.layers = 30;
box = [1.564 1.21 1.1235];
opt.box=box(1:2);

x = [box(1)*rand(N,1) box(2)*rand(N,1) box(3)*rand(N,1)];
q = rand(N,1); q = q-mean(q);
idx = 1;

xi = linspace(2,10,5);

for j=1:length(xi)
    opt.xi = xi(j);
    phi1 = SE2P_direct_real_mex(idx,x,q,opt);
    phi2 = SE2P_direct_fd_mex(idx,x,q,opt);
    phi3 = SE2P_direct_k0_mex(idx,x,q,opt);
    phi4 = SE2P_direct_self_mex(idx,x,q,opt);
    phi(j) = phi1+phi2+phi3+phi4;
end
phi(:)

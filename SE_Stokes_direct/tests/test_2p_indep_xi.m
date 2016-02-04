clear all

N = 10;
opt.layers = 30;
box = [1.564 1.21 1.1235];
opt.box=box(1:2);

x = [box(1)*rand(N,1) box(2)*rand(N,1) box(3)*rand(N,1)];
f = rand(N,3); f = f-repmat(mean(f),N,1);
idx = 1;

xi = linspace(2,10,5);

for j=1:length(xi)
    opt.xi = xi(j);
    u1 = SE2P_Stokes_direct_real_mex(idx,x,f,opt);
    u2 = SE2P_Stokes_direct_fd_mex(idx,x,f,opt);
    u3 = SE2P_Stokes_direct_k0_mex(idx,x,f,opt);
    u4 = SE2P_Stokes_direct_self_mex(idx,x,f,opt);
    u(j,:) = u1+u2+u3+u4;
end

e = diff(u, 1, 1);
e_max = norm(e(:), inf);
assert(e_max < 1e-12, '2P Stokes xi independence failed')
disp('PASSED')
clear
rng(1)
N = 2;
L = 1;
box = [L L L];
[x, f] = vector_system(N, box);
opt.box = box;
opt.layers = 100;

xi = linspace(1,5,5);
idx = 1;

for j=1:length(xi)
    opt.xi = xi(j);
    opt.rc = 1;
    phi1 = SE2P_direct_real_mex(idx,x,f,opt);
    phi2 = SE2P_direct_fd_mex(idx,x,f,opt);
    phi3 = SE2P_direct_k0_mex(idx,x,f,opt);
    phi4 = SE2P_direct_self_mex(idx,x,f,opt);
    phi(j) = phi1+phi2+phi3+phi4;
end

err = phi(1)-phi(2:end);
rms_err = rms(err) / rms(phi(1))

if rms_err < 1e-14
    fprintf('\n********** XI INDEPENDENCE: OK **********\n\n')
else
    error('XI INDEPENDENCE: FAILED')
end

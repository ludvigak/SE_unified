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

MATLAB = false;

F = zeros(length(xi),3);
for j=1:length(xi)
    opt.xi = xi(j);
    opt.rc = 6/opt.xi;
    F1 = se1p_real_space_force(idx,x,f,opt,MATLAB);
    F2 = se1p_fourier_space_direct_force(idx,x,f,opt,MATLAB);
    F3 = se1p_k0_direct_force(idx,x,f,opt,MATLAB)

    F(j,:) = F1+F2+F3;
end

err = bsxfun(@minus,F(1,:),F(2:end,:));

rms_err = rms(err(:)) / rms(F(1,:))

if rms_err < 1e-14
    fprintf('\n********** XI INDEPENDENCE: OK **********\n\n')
else
    error('XI INDEPENDENCE: FAILED')
end
    
clear
% Test Laplace Spectral Ewald
rng(1);
N = 200;
L = 3;
box = [L L L];
[x, f] = laplace_system(N, box);


% Try integer and non-integer oversampling
warning('OFF', 'FSE:OversamplingIncreased')
for s = [3 (1+sqrt(3))] 
    % Try odd and even base grid
    for M0 = [10 11]
        opt.M = M0*box;
        opt.xi = pi*M0 / 12;
        opt.P = 24;
        opt.oversampling = s;
        opt.rc = 6 / opt.xi;
        opt.box = box;

        % Ewald
        pre = laplace_precomp(opt);
        uf = laplace_fourier_space(x, f, opt, pre);
        ur = laplace_real_space(x, f, opt);
        us = -f*opt.xi *2/sqrt(pi);
        ue = uf+ur+us;
        % Direct
        u = laplace_direct(x, f, box);

        error = norm(u-ue, inf) / norm(u,inf);
        assert(error < 1e-13, 'Laplace SE did not match direct sum');
    end
end

function G = laplace_k_scaling(H, fse_opt, opt, pre)

G = cell(1);

% k-vectors
[k1 k2 k3] = k_vectors(fse_opt.padded_M, fse_opt.padded_box);
[K1 K2 K3] = ndgrid(k1,k2,k3);
Ksq = K1.^2 + K2.^2 + K3.^2;

% scale
K = sqrt(Ksq);
B = exp(-(1-opt.eta)*Ksq/(4*opt.xi^2));
G{1} = H{1} .* B .* pre.GR;
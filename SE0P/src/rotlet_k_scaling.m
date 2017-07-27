function G = rotlet_k_scaling(H, fse_opt, opt, pre)

G = cell(3, 1);
GR = -pre.GR; % why minus here? Because of the code!


% Pure MATLAB:
% k-vectors
[k1 k2 k3] = k_vectors(fse_opt.padded_M, fse_opt.padded_box);
[K1 K2 K3] = ndgrid(k1,k2,k3);
Ksq = K1.^2 + K2.^2 + K3.^2;

% scale
B = exp(-(1-opt.eta)*Ksq/(4*opt.xi^2)) .* pre.GR / (4*pi); % keep it real
B(k1==0, k2==0, k3==0) = 0;
% G = B*(HxK)
G{1} = -1i*4*pi*B.*(H{2}.*K3 - H{3}.*K2);
G{2} = -1i*4*pi*B.*(H{3}.*K1 - H{1}.*K3);
G{3} = -1i*4*pi*B.*(H{1}.*K2 - H{2}.*K1);



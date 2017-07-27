function G = stokeslet_k_scaling(H, fse_opt, opt, pre)

G = cell(3, 1);
GR = -pre.GR; % why minus here? Because the paper says so!

% MEX
[G{1} G{2} G{3}] = SE0P_stokeslet_fast_fs_k_scaling(H{1}, H{2}, H{3}, GR, opt.xi, ...
						   fse_opt.padded_box, opt.eta);

% % Pure MATLAB:
% % k-vectors
% [k1 k2 k3] = k_vectors(fse_opt.padded_M, fse_opt.padded_box);
% [K1 K2 K3] = ndgrid(k1,k2,k3);
% Ksq = K1.^2 + K2.^2 + K3.^2;

% % scale
% B = (1 + Ksq/(4*opt.xi^2)) .* GR .* exp(-(1-opt.eta)*Ksq/(4*opt.xi^2));
% KdotH = K1.*H{1} + K2.*H{2} + K3.*H{3};
% G{1} = B .* (Ksq.*H{1} - KdotH.*K1);
% G{2} = B .* (Ksq.*H{2} - KdotH.*K2);
% G{3} = B .* (Ksq.*H{3} - KdotH.*K3);



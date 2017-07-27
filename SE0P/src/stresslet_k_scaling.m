function G = stresslet_k_scaling(H, fse_opt, opt, pre)

G = cell(3, 1);
GR = pre.GR; 

% - Scaling is missing a minus sign
% - Code needs a minus sign
% so this is right

% MEX scaling
[G{1} G{2} G{3}] = SE0P_stresslet_fast_fs_k_scaling(H, GR, opt.xi, ...
						   fse_opt.padded_box, opt.eta);
return

% Pure MATLAB:
% k-vectors
[k1 k2 k3] = k_vectors(fse_opt.padded_M, fse_opt.padded_box);
K = cell(3,1);
[K{1} K{2} K{3}] = ndgrid(k1,k2,k3);
Ksq = K{1}.^2 + K{2}.^2 + K{3}.^2;
xi = opt.xi;
% Beenakker
%B = -1i*GR .* ...
%    (1 + Ksq/(4*xi^2) + Ksq.^2/(8*xi^4)) ...
%    .* exp(-(1-opt.eta)*Ksq/(4*xi^2));
% Hasimoto
B = -1i*GR .* ...
    (1 + Ksq/(4*xi^2)) ...
    .* exp(-(1-opt.eta)*Ksq/(4*xi^2));
for j=1:3
    tmp = zeros(fse_opt.padded_M);
    for l=1:3
        for m=1:3
            tmp = tmp + ...
                  (Ksq.*( (j==l)*K{m} + (l==m)*K{j} + (m==j)*K{l} ) ...
                   - 2*K{j}.*K{l}.*K{m}) ...
                  .* H{l,m};            
            
        end
    end
    G{j} = B.*tmp;
end

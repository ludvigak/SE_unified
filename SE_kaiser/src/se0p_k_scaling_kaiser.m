function G = se0p_k_scaling_kaiser(H, fse_opt, pre_window, pre_kernel)

% k-vectors
[k1 k2 k3] = k_vectors(fse_opt.padded_M, fse_opt.padded_box);
[K1 K2 K3] = ndgrid(k1,k2,k3);
Ksq = K1.^2 + K2.^2 + K3.^2;

% scale
B = exp(-Ksq/(4*fse_opt.xi^2));
G = H .* B .* pre_kernel.GR .* pre_window.F;


% ----------------------------------------
function [k1 k2 k3] = k_vectors(M,box)
k1 = k_vec(M(1), box(1));
k2 = k_vec(M(2), box(2));
k3 = k_vec(M(3), box(3));

% ----------------------------------------
function [k] = k_vec(M,L)
if mod(M,2)==0
    MM = M/2;
    k = (2*pi/L)*[-MM:(MM-1)];
elseif mod(M-1,2)==0
    MM = (M-1)/2;
    k = (2*pi/L)*[-MM:MM];
else error('k-vectors not computed');
end

% k = fftshift(k);
% idx = 1;
% if(abs(k(idx))>eps)
%     k = circshift(k,1);
% end

% assert(abs(k(idx))<eps)



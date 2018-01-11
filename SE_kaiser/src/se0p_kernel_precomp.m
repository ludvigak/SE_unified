function pre = se0p_kernel_precomp(fse_opt)
%% Pecomputation for resolving the free-space kernel 

% Modified Greens function, corresponding to delta in corner
[k1 k2 k3] = k_vectors(fse_opt.oversampled_M, fse_opt.oversampled_box);
[K1 K2 K3] = ndgrid(k1,k2,k3);
K = sqrt(K1.^2 + K2.^2 + K3.^2);
GR = 8*pi*(sin(K*fse_opt.R/2)./K).^2;
GR(K==0) = 2*pi*fse_opt.R^2;
GR = ifftn(GR);
% GR now has rubbish in the center,
% truncate in real space by shifting, picking out center and shifting back.
% This is faster than centering by multiplying with e^(i*k*xc) in k-space.
n = fse_opt.padded_M;
N = fse_opt.oversampled_M;
GR = fftshift(GR);
start = floor( (N-n)/2 );
GR = GR( start(1) + (1:n(1)), ...
         start(2) + (1:n(2)), ...
         start(3) + (1:n(3)) );
GR = ifftshift(GR);
pre.GR = fftn(GR);


% ---------------------------------------------------------
function [k1 k2 k3] = k_vectors(M,box)
k1 = k_vec(M(1), box(1));
k2 = k_vec(M(2), box(2));
k3 = k_vec(M(3), box(3));

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
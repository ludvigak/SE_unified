function G = se3p_k_scaling_kaiser(G,pre,opt)
    
eta = opt.eta; xi = opt.xi;

% k-vectors
[k1 k2 k3] = k_vectors(opt.M,opt.box); 
[K1 K2 K3] = ndgrid(k1,k2,k3);

% scale
k2 = K1.^2 + K2.^2 + K3.^2;

Z = exp(-k2/(4*xi^2))./k2.*pre.F;
Z(1,1,1) = 0;
G = G.*Z;

% ===========================================
function [k1 k2 k3] = k_vectors(M,box)

K=[];
for i=1:3
    if mod(M(i),2)==0
        MM = M(i)/2;
        k = (2*pi/box(i))*[0:(MM-1) -MM:-1];
    else
        MM = (M-1)/2;
        k = (2*pi/box(i))*[0:MM -MM:-1];
    end
    K = [K; k];
end
k1 = K(1,:);
k2 = K(2,:);
k3 = K(3,:);

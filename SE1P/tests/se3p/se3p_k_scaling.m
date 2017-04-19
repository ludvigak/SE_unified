function G = se3p_k_scaling(G, opt)

eta = opt.eta;
xi = opt.xi;
M = opt.M;

% k-vectors
[k1, k2, k3] = k_vec(M,opt.box);
[K1, K2, K3] = ndgrid(k1,k2,k3);

% scale
k2 =  K1.^2 + K2.^2 + K3.^2;
Z = exp(-(1-eta)/(4*xi^2)*k2)./k2;
Z(1,1,1) = 0;

G = G.*Z;

end

% ------------------------------------------------------------------------------
function [k1, k2, k3] = k_vec(M,box)
    if (all(mod(M,2)==0))
        MM = M/2;
        k1 = (2*pi/box(1))*[0:(MM(1)-1) -MM(1):-1];
        k2 = (2*pi/box(2))*[0:(MM(2)-1) -MM(2):-1];
        k3 = (2*pi/box(3))*[0:(MM(3)-1) -MM(3):-1];
        
    elseif(all(mod(M-1,2)==0))
        MM = (M-1)/2;
        k1 = (2*pi/box(1))*[0:MM(1) -MM(1):-1];
        k2 = (2*pi/box(2))*[0:MM(2) -MM(2):-1];
        k3 = (2*pi/box(3))*[0:MM(3) -MM(3):-1];
        
    else error('k-vectors not computed (FIXME)');
    end
end
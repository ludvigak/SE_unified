function [G, Gr, G0] = se2p_k_scaling(G,Gr,G0,opt)

R = opt.R; eta = opt.eta; xi = opt.xi;

% k-vectors
[k1 zidx1] = k_vectors(opt.M, opt.L, 1);
[k2 zidx2] = k_vectors(opt.M, opt.L, 1);
kappa      = k_vectors(opt.Mz,opt.Lz,1);
[K1 K2 KAPPA] = ndgrid(k1,k2,kappa);

% scale the whole domain
ksq = K1.^2 + K2.^2 + KAPPA.^2;
Znum = exp(-(1-eta)/(4*xi^2)*ksq);
Z = Znum./ksq;
Z(zidx1,zidx2,1) = 0;

% scale zero mode
kappa = k_vectors(opt.Mz, opt.Lz, opt.s0)';
ksq = k1(zidx1)^2 + k2(zidx2)^2 + kappa.^2;
kmod  = sqrt(ksq);
Znum = exp(-(1-eta)/(4*xi^2)*ksq);

Green=-1./kmod.^2.*(R*kmod.*sin(R*kmod)+cos(R*kmod)-1);
Z0 = Znum.*Green;

% % Finite limit at k3=0.
Z0(1) = R^2/2;

% scale the local pad
if(numel(opt.local_pad)>0)
    % local pad scaling
    kappa   = k_vectors(opt.Mz, opt.Lz, opt.s);
    
    [K1,K2, KAPPA] = ndgrid(k1(opt.local_pad),k2(opt.local_pad),kappa);
    ksq = K1.^2 + K2.^2 + KAPPA.^2;
    Znum = exp(-(1-eta)/(4*xi^2)*ksq);

    Zr = Znum./ksq;
    Zr(isinf(Zr)==1)=0;
end

G0 = Z0.*G0;
Gr = Zr.*Gr;
G  = Z.*G;

end

% ------------------------------------------------------------------------------
function [ks, idxz] = k_vectors(M,L,q)

% q: oversampling factor (1 = no oversampling)
Mq = round(q*M);

if mod(Mq,2)==0
    k = (-Mq/2):(Mq/2-1);
else
    k = -(Mq-1)/2:(Mq-1)/2;
end

k = fftshift(k); % standard reordering

idxz = 1;
ks = 2*pi*k/(L*q);
if (abs(ks(idxz))>eps)
    ks = circshift(ks,[1 1]);
end

assert(abs(ks(idxz))<eps)
end

function [G, Gres, G0] = se1p_k_scaling(G,Gres,G0,opt)

R = opt.R; eta = opt.eta; xi = opt.xi;

% k-vectors
kappa_1   = k_vectors(opt.Mx, opt.Lx, opt.sg);
kappa_2   = k_vectors(opt.My, opt.Ly, opt.sg);
[k, zidx] = k_vectors(opt.M,  opt.L,  1); % z direction is periodic;no oversampling
[KAPPA1, KAPPA2, K] = ndgrid(kappa_1,kappa_2,k);

% scale the whole domain
k2 =  KAPPA1.^2 + KAPPA2.^2 + K.^2;
Znum = exp(-(1-eta)/(4*xi^2)*k2);
Z = Znum./k2;
Z(1,1,1) = 0;

if(opt.sg~=opt.s0)
    kappa_1   = k_vectors(opt.Mx, opt.Lx, opt.s0);
    kappa_2   = k_vectors(opt.My, opt.Ly, opt.s0);
    [KAPPA1, KAPPA2, K] = ndgrid(kappa_1,kappa_2,k);
    k2 =  KAPPA1.^2 + KAPPA2.^2 + K.^2;
    Znum = exp(-(1-eta)/(4*xi^2)*k2);
end

kmod  = sqrt(k2(:,:,opt.k0mod));
Green=(1-besselj(0,R*kmod))./k2(:,:,zidx)-R*log(R)*besselj(1,R*kmod)./kmod;
Z0 = Znum(:,:,opt.k0mod).*Green;
% Finite limit at k3=0.
Z0(1,1) = R^2/4-R*log(R)*R/2;

if(opt.sg==opt.s0 && opt.sl==opt.s0)
    Z(:,:,opt.k0mod) = Z0; 
    G = Z.*G;
    return
elseif(opt.sg~=opt.s0)
    G0 = Z0.*G0;
    G = Z.*G;
end

if(opt.sg~=opt.sl && numel(opt.local_pad)>0)
    % local pad scaling
    kappa_1   = k_vectors(opt.Mx, opt.Lx, opt.sl);
    kappa_2   = k_vectors(opt.My, opt.Ly, opt.sl);
    
    [KAPPA1, KAPPA2, K] = ndgrid(kappa_1,kappa_2,k(opt.local_pad));
    k2 = KAPPA1.^2 + KAPPA2.^2 + K.^2;
    Znum = exp(-(1-eta)/(4*xi^2)*k2);
    Zres = Znum./k2;
    Gres = Gres.*Zres;
end


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
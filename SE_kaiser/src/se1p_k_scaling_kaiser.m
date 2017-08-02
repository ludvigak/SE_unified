function [G, Gr, G0] = se1p_k_scaling_kaiser(G,Gr,G0,pre,opt)

R = opt.R; xi = opt.xi;

% k-vectors
[k1 zidx] = k_vectors(opt.M, opt.L, 1);
kappa_1   = k_vectors(opt.My,opt.Ly,1);
kappa_2   = k_vectors(opt.Mz,opt.Lz,1);
[K1 KAPPA1 KAPPA2] = ndgrid(k1,kappa_1,kappa_2);

% scale the whole domain
ksq = K1.^2 + KAPPA1.^2 + KAPPA2.^2;
Znum = exp(-ksq/(4*xi^2));
Z = Znum./ksq.*pre.F;
Z(zidx,1,1) = 0;

% scale zero mode
kappa_1 = k_vectors(opt.My, opt.Ly, opt.s0);
kappa_2 = k_vectors(opt.Mz, opt.Lz, opt.s0);
[KAPPA1 KAPPA2] = meshgrid(kappa_1,kappa_2);
ksq = k1(zidx)^2+KAPPA1.^2+KAPPA2.^2;
kmod= sqrt(ksq);

Znum = exp(-ksq/(4*xi^2));
Green=(1-besselj(0,R*kmod))./ksq-R*log(R)*besselj(1,R*kmod)./kmod;

% Finite limit at k1=0.
Z0 = Znum.*Green.*pre.F0;
Z0(1,1)=(R^2/4-R*log(R)*R/2);


% scale the local pad
if(numel(opt.local_pad)>0)
    % local pad scaling
    kappa_1   = k_vectors(opt.My, opt.Ly, opt.s);
    kappa_2   = k_vectors(opt.Mz, opt.Lz, opt.s);
    
    [K1,KAPPA1, KAPPA2] = ndgrid(k1(opt.local_pad),kappa_1,kappa_2);
    ksq = K1.^2 + KAPPA1.^2 + KAPPA2.^2;
    Zr = exp(-ksq/(4*xi^2)).*pre.Fr./ksq;
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

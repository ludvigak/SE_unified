function [G, Gres, G0, varargout] = se1p_k_scaling(G,Gres,G0,opt)

R = opt.R; eta = opt.eta; xi = opt.xi;

% k-vectors
[k, zidx] = k_vectors(opt.M,  opt.L,  1); % z direction is periodic;no oversampling
kappa_1   = k_vectors(opt.My, opt.Ly, opt.sg);
kappa_2   = k_vectors(opt.Mz, opt.Lz, opt.sg);
[K, KAPPA1, KAPPA2] = ndgrid(k,kappa_1,kappa_2);

% scale the whole domain
ksq =  K.^2 + KAPPA1.^2 + KAPPA2.^2;
Znum = exp(-(1-eta)/(4*xi^2)*ksq);
Z = Znum./ksq;
Z(zidx,1,1) = 0;

scale_t = tic;
if(opt.sg~=opt.s0)
    kappa_1   = k_vectors(opt.My, opt.Ly, opt.s0);
    kappa_2   = k_vectors(opt.Mz, opt.Lz, opt.s0);
    [K, KAPPA1, KAPPA2] = ndgrid(k,kappa_1,kappa_2);
    ksq =  K.^2 + KAPPA1.^2 + KAPPA2.^2;
    ksq = squeeze(ksq(opt.k0mod,:,:));
    Znum = exp(-(1-eta)/(4*xi^2)*ksq);
end
walltime = toc(scale_t);

kmod  = sqrt(ksq);
Green=(1-besselj(0,R*kmod))./ksq-R*log(R)*besselj(1,R*kmod)./kmod;
Z0 = Znum.*Green;

% Finite limit at k3=0.
Z0(1,1) = R^2/4-R*log(R)*R/2;
whos
if(opt.sg==opt.s0 && opt.sl==opt.s0)
    Z(opt.k0mod,:,:) = Z0; 
    G = Z.*G;
    return
elseif(opt.sg~=opt.s0)
    G0 = Z0.*G0;
    G = Z.*G;
end
   
if(opt.sg~=opt.sl && numel(opt.local_pad)>0)
    % local pad scaling
    kappa_1   = k_vectors(opt.My, opt.Ly, opt.sl);
    kappa_2   = k_vectors(opt.Mz, opt.Lz, opt.sl);
    
    [K, KAPPA1, KAPPA2] = ndgrid(k(opt.local_pad),kappa_1,kappa_2);
    ksq = K.^2 + KAPPA1.^2 + KAPPA2.^2;
    Znum = exp(-(1-eta)/(4*xi^2)*ksq);

    scale_t = tic;
    Zres = Znum./ksq; 
    Zres(isinf(Zres)==1)=0;
    Gres = Gres.*Zres;
    walltime = walltime + toc(scale_t);
end

if nargout == 4
    varargout{1} = walltime;
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

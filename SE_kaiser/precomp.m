function pre = precomp(opt)

    L = opt.box(1);
    M = opt.M(1);
    h = L/M;
    P = opt.P;
    w = P/2;
    beta   = opt.beta;

    
    x = 0:h:L-h;
    f = kaiser((x-L/2)/h,beta,w)*h;
    
    F=fft(f);
    [F1 F2 F3] = ndgrid(F.^2);
    F = F1.*F2.*F3;
    pre.F = 1./F;
    pre.F(isnan(F))=0;

function z= kaiser(x,beta,w)

    t = sqrt(1.-(x/w).^2);
    z = exp(beta*(t-1)).*(abs(x)<=w);
    
    
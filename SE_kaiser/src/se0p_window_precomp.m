function pre = se0p_window_precomp(opt)

    L=opt.oversampled_box(1);
    M=opt.oversampled_M(1);
    h = L/M;
    P = opt.P;
    w = P/2;
    s = opt.oversampling;s=1;
    beta = opt.beta;

    % common in all parts Lx=L
    x = 0:h:s*L-h;
    f = kaiser((x-s*L/2)/h,beta,w)*h;
    n = opt.padded_M(1);
    N = opt.oversampled_M(1);
    start = round( (N-n)/2 );
    f = f(start + (1:n));
    F = fft(f);
    F2= 1./F.^2;
    [F1 F2 F3] = ndgrid(F2, F2, F2);
    F = F1.*F2.*F3;
    pre.F = F;

function z= kaiser(x,beta,w)

    t = sqrt(1.-(x/w).^2);
    z = exp(beta*(t-1)).*(abs(x)<=w);
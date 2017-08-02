function pre = se1p_window_precomp(opt)
  
    % We assume periodicity in x and y and free in z.
    
    L = opt.box(1);
    Ly= opt.Ly;
    Lz= opt.Lz;
    M = opt.M;
    Mz= opt.Mz;
    My= opt.My;
    h = L/M;
    P = opt.P;
    w = P/2;
    s = opt.s;
    n = opt.local_pad;
    s0= opt.s0;
    beta = opt.beta;

    % common in all parts Lx=L
    x = 0:h:L-h;
    f = kaiser((x-L/2)/h,beta,w)*h;
    F=fft(f);          F2=F.^2;
    
    % global set with no upsampling
    z = 0:h:Lz-h;
    f = kaiser((z-Lz/2)/h,beta,w)*h;
    Fz=fft(f,Mz);  Fz2=Fz.^2;
    y = 0:h:Ly-h;
    g = kaiser((y-Ly/2)/h,beta,w)*h;
    Fy=fft(g,My);  Fy2=Fy.^2;
    [f1 f2 f3] = ndgrid(F2,Fy2,Fz2);
    F = f1.*f2.*f3;
    pre.F = 1./F;
%    pre.F(:,My/2+1,Mz/2+1) = 0;
    
    % local set with s upsampling on local pad
    z = 0:h:s*Lz-h;
    f = kaiser((z-s*Lz/2)/h,beta,w)*h;
    Fz=fft(f,round(s*Mz));  Fz2=Fz.^2;
    y = 0:h:s*Ly-h;
    g = kaiser((y-s*Ly/2)/h,beta,w)*h;
    Fy=fft(g,round(s*My));  Fy2=Fy.^2;
    [f1 f2 f3] = ndgrid(F2(n),Fy2,Fz2);
    Fr = f1.*f2.*f3;
    pre.Fr = 1./Fr;
    
    % zero mode with s0 upsampling
    z = 0:h:s0*Lz-h;
    f = kaiser((z-s0*Lz/2)/h,beta,w)*h;
    Fz=fft(f,round(s0*Mz)); Fz2=Fz.^2;
    y = 0:h:s0*Ly-h;
    g = kaiser((y-s0*Ly/2)/h,beta,w)*h;
    Fy=fft(g,round(s0*My)); Fy2=Fy.^2;
    [f2 f3] = meshgrid(Fy2,Fz2);
    F0 = F2(1)*f2.*f3;
    pre.F0 = 1./F0;

function z= kaiser(x,beta,w)

    t = sqrt(1.-(x/w).^2);
    z = exp(beta*(t-1)).*(abs(x)<=w);
    
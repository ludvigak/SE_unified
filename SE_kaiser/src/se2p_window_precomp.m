function pre = se2p_window_precomp(opt)
   
    % We assume periodicity in x and y and free in z.
    
    L = opt.box(1);
    Lz= opt.Lz;
    M = opt.M;
    Mz= opt.Mz;
    h = L/M;
    P = opt.P;
    w = P/2;
    s = opt.s;
    n = opt.local_pad;
    s0= opt.s0;
    beta = opt.beta;
   
    % common in all parts Lx=Ly=L
    x = 0:h:L-h;
    f = kaiser((x-L/2)/h,beta,w)*h;
    F=fft(f);          F2=F.^2;
    
    % global set with no upsampling
    x = 0:h:Lz-h;
    f = kaiser((x-Lz/2)/h,beta,w)*h;
    Fz=fft(f,Mz);  Fz2=Fz.^2;
    [f1 f2 f3] = ndgrid(F2,F2,Fz2);
    F = f1.*f2.*f3;
    pre.F = 1./F;
    pre.F(:,:,Mz/2+1) = 0;
    
    % local set with s upsampling on local pad
    x = 0:h:s*Lz-h;
    f = kaiser((x-s*Lz/2)/h,beta,w)*h;
    Fz=fft(f,round(s*Mz));  Fz2=Fz.^2;
    [f1 f2 f3] = ndgrid(F2(n),F2(n),Fz2);
    Fr = f1.*f2.*f3;
    pre.Fr = 1./Fr;
    %    pre.Fr(:,:,round(s*Mz)/2+1) = 0;
    
    
    % zero mode with s0 upsampling
    x = 0:h:s0*Lz-h;
    f = kaiser((x-s0*Lz/2)/h,beta,w)*h;
    Fz=fft(f',round(s0*Mz));
    Fz2=Fz.^2;
    f1 = F2(1);f2=F2(1);f3=Fz2;
    F0 = f1*f2*f3;
    pre.F0 = 1./F0;

%pre.F0(round(s0*Mz/2)+1)=0;

function z= kaiser(x,beta,w)

    t = sqrt(1.-(x/w).^2);
    z = exp(beta*(t-1)).*(abs(x)<=w);
    
    
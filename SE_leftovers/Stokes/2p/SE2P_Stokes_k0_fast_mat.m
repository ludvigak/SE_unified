function G = SE2P_Stokes_k0_fast_mat(idx,x,opt)
    
    verb = true;
    
    % unpack parameters
    M = opt.k0_M;
    xi = opt.xi;
    box = opt.box;
    
    % system size
    x = x(idx,:);
    N = size(x,1);
    
    cprintf(verb,'[SE2P Stokes (K0F)] Matrix form. N=%d, M=%d... ',N,M);    
    assert(N*M<1e8,'System too large (%d,%d)',N,M)
    
    % gauss/lobotto points
    z_gauss = gauss_points(0,box(3),M); 
    z_x = x(:,3);
    
    G = zeros(M,N);
    for k=1:M
        znk=z_x-z_gauss(k);
        G(k,:)=exp(-xi^2*znk.^2)/(2*xi) + sqrt(pi)*znk.*erf(xi*znk);
    end
    
    cprintf(verb,'Done\n');
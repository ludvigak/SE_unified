clear all, close all

s='matform';

switch s
  case 'err'
    N = 1000;
    box = [1 1 1]
    [x f] = SE2P_state(N,box);
    opt.box = box(1:2);

    opt.xi=10;
    
    uref = SE2P_Stokes_direct_k0_mex(1:N,x,f,opt);
    
    opt.box = box;
    opt.k0_M=100;
    
    u=SE2P_Stokes_k0_fast(1:N,x,f,opt);
    
    max(abs(u(:)-uref(:)))    
    sqrt((1/N)*sum( (u-uref).^2,1) )
    
  case 'perf'
    
    N=500000;
    box = [1 1 1]
    [x f] = SE2P_state(N,box);
    opt.box = box(1:2);
    opt.xi=10;
    opt.box = box;
    opt.k0_M=150;
    
    tic
    u=SE2P_Stokes_k0_fast(1:N,x,f,opt);
    toc
    
  case 'matform'
    
    N=500000;
    box = [1 1 1];
    [x f] = SE2P_state(N,box);
    opt.box = box(1:2);
    opt.xi=10;
    opt.box = box;
    opt.k0_M=150;
    
    tic
    G=SE2P_Stokes_k0_fast_mat(1:N,x,opt);
    toc

    tic
    u=SE2P_Stokes_k0_fast_apply(1:N,x,f,G,opt);
    u=SE2P_Stokes_k0_fast_apply(1:N,x,f,G,opt);
    u=SE2P_Stokes_k0_fast_apply(1:N,x,f,G,opt);
    toc

  case 'matform err'
        
    N=1000;
    box = [1 1 1];
    [x f] = SE2P_state(N,box);
    opt.box = box(1:2);
    opt.xi=10;
    opt.box = box;
    opt.k0_M=100;

    uref = SE2P_Stokes_direct_k0_mex(1:N,x,f,opt);
    
    G=SE2P_Stokes_k0_fast_mat(1:N,x,opt);
    u=SE2P_Stokes_k0_fast_apply(1:N,x,f,G,opt);

    max(abs(u(:)-uref(:)))    
    sqrt((1/N)*sum( (u-uref).^2,1) )

end
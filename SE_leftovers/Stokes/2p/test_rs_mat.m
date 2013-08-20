clear all, close all

s = 'rc sph periodic'

switch s 
  case 'full'

    N = 2000;
    box = [1 1 1];
    xi = 2;
    
    [x f] = SE2P_state(N,box);
    F = f(:);
    
    A = SE_Stokes_rs_mat(x,-1,xi);
    U = A*F;

    % reference computation
    opt.xi = xi;
    opt.box=box(1:2);
    opt.layers=0; % corresponds to matrix form only for zero layers
    u=SE2P_Stokes_direct_real_mex(1:N,x,f,opt);

    max(abs(U-u(:)))
    
  case 'rc ref'
    
    N = 1000;
    box = [1 1 1];
    xi = 2;
    
    [x f] = SE2P_state(N,box);
    F = f(:);
    
    rc = 1/5;
    
    A = SE_Stokes_rs_mat(x,rc,xi);
    U = A*F;

    % reference computation
    opt.xi = xi;
    opt.box=box(1:2);
    opt.layers=0; % corresponds to matrix form only for zero layers
    opt.rc = rc;
    u=SE2P_Stokes_direct_rsrc_mex(1:N,x,f,opt);

    max(abs(U-u(:)))
    
  case 'rc prof'
    
    N = 20000;
    box = [1 1 1];
    xi = 2;
    
    [x f] = SE2P_state(N,box);
    F = f(:);
    
    rc = 1/10;    
    A = SE_Stokes_rs_mat(x,rc,xi);
    U = A*F;    
    
  case 'rc sph periodic'
    xi = 2;
    box = [1 1 1];
    N_sph = 500;
    r_sph=1/30;
    d_min = 3*r_sph;
    rc = 1/5;
    
    % place spheres with specified minimal distance.
    % also, find neghbour spheres (counting 2-periodicity)
    [cen sph_nbh nbh_shift] = sph_coords(N_sph, r_sph, d_min, box(1),rc);

    % collect points 
    x = [];
    for k=1:N_sph
        xx = sph_pts(r_sph,cen(k,:),'icosa',0);
        x = [x; xx];
    end
    M = size(xx,1); % points per sphere
    N = size(x,1);
    f = rand(size(x));
    
    % reference computation: one layer, enforce rc
    opt.xi = xi;
    opt.box=box(1:2);
    opt.layers=1; 
    opt.rc = rc;
    u=SE2P_Stokes_direct_rsrc_mex(1:N,x,f,opt);

    % matrix form, assembly guided by neghbour spheres.  
    A = SE_Stokes_rs_mat(x, rc, opt.xi, sph_nbh, M, nbh_shift);

    % apply
    U = A*f(:);
    
    % error (relative)
    max(abs(  (U-u(:))./U ) )
    
end

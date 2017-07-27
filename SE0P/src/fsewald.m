function fse = fsewald(kernel)

switch kernel
  case 'laplace'
    generator = @laplace_system;
    direct_sum = @laplace_direct;
    self = @(f, opt) -f*opt.xi *2/sqrt(pi);
    real_sum = @laplace_real_space;    
    fourier_sum = @laplace_fourier_space;
    fourier_sum_3p = @spectral_ewald;   
    fourier_estimate = @estimates.laplace.fourier;
    real_estimate = @() inf;;    
    precomp = @laplace_precomp;
    calc_M  = @estimates.laplace.calc_M;
    calc_P  = @estimates.laplace.calc_P;
    calc_RC = @estimates.laplace.calc_rc;
  case 'stokeslet'
    generator = @vector_system;
    direct_sum = @stokeslet_direct;
    self = @(f, opt) -4*opt.xi/sqrt(pi)*f;
    real_sum = @stokeslet_real_space;
    fourier_sum = @stokeslet_fourier_space;
    fourier_sum_3p = @SE_Stokes;
    fourier_estimate = @estimates.stokeslet.fourier;
    real_estimate = @estimates.stokeslet.real;
    precomp = @stokeslet_precomp;
    calc_M  = @estimates.stokeslet.calc_M;
    calc_P  = @estimates.laplace.calc_P;
    calc_RC = @estimates.stokeslet.calc_rc;
  case 'rotlet'
    generator = @vector_system;
    direct_sum = @rotlet_direct;
    self = @(f, opt) 0;
    real_sum = @rotlet_real_space;
    fourier_sum = @rotlet_fourier_space;
    fourier_sum_3p = @(idx,x,f,xi,opt) SE_Rotlet(x(idx,:),x,f,xi,opt);
    fourier_estimate = @estimates.rotlet.fourier;
    real_estimate = @estimates.rotlet.real;
    precomp = @rotlet_precomp;     
    calc_M  = @estimates.rotlet.calc_M;
    calc_P  = @estimates.rotlet.calc_P;
    calc_RC = @estimates.rotlet.calc_rc; 
  case 'stresslet'
    generator = @(N, box) vector_system(N, box, 6);
    direct_sum = @stresslet_direct;
    self = @(f, opt) 0;
    real_sum = @stresslet_real_space;
    fourier_sum = @stresslet_fourier_space;
    fourier_sum_3p = @SE_Stresslet;
    fourier_estimate = @estimates.stresslet.fourier;
    real_estimate = @estimates.stresslet.real;    
    precomp = @stresslet_precomp;     
    calc_M  = @estimates.stresslet.calc_M;
    calc_P  = @estimates.laplace.calc_P;
    calc_RC = @estimates.stresslet.calc_rc;   
  otherwise
    error('Unknown kernel')
end

fse.generator       = generator      ;
fse.direct_sum      = direct_sum     ;
fse.self            = self           ;
fse.real_sum        = real_sum       ;
fse.fourier_sum     = fourier_sum    ;
fse.fourier_sum_3p  = fourier_sum_3p ;
fse.fourier_estimate= fourier_estimate;
fse.real_estimate   = real_estimate  ;
fse.precomp         = precomp        ;
fse.calc_M          = calc_M         ;
fse.calc_P          = calc_P         ;
fse.calc_RC         = calc_RC        ;

 









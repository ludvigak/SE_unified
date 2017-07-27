clear

tol = 10.^(0:-1:-13);
kernels = {'stokeslet', 'stresslet', 'rotlet'};


col = {'b','r','k'};
rep = 4;

N = 20000;
L = 2;

box = [L L L];
opt.box = box;
opt.xi = 7;
opt.oversampling = 1+sqrt(3);
opt.no_extra_support = true;

rms = @(u) sqrt(1/size(u,1)*sum(u(:).^2));


for idx_kernel = 1:numel(kernels)
    kernel = kernels{idx_kernel};
    fse = fsewald(kernel);
    rng(1);
    [x f] = fse.generator(N, box);
    
    if strcmp(kernel,'stokeslet')
        fse.real_sum = @stokeslet_real_space_rangesearch;
    end
    switch kernel
        % Optimized at N = 20000 for error < .51e-8
      case 'stokeslet'
        EM = -4;
        EP = -2;
        Erc = -.1;
      case 'stresslet'
        EP = -2;
        EM = -6;
        Erc = -.15;
      case 'rotlet'
        EM = 0;
        EP = -2;
        Erc = -.1;
    end
    
    errors = [];
    uref = fse.direct_sum(x, f, box);
    
    for idx_tol = 1:numel(tol);
        Tol = tol(idx_tol);

        opt.M = fse.calc_M(f,opt,Tol)*[1 1 1] + EM;
	if(idx_tol==12 && ( strcmp(kernel,'stokeslet') || strcmp(kernel,'stresslet') ) )
		opt.M = opt.M + 2;
	end
        opt.P = fse.calc_P(f,opt,Tol) + EP;
        opt.rc= fse.calc_RC(f,opt,Tol) + Erc;
        
        times.pre(idx_tol) = inf;
        times.real(idx_tol) = inf;
        times.fourier(idx_tol) = inf;
        times.fft(idx_tol) = inf;
        times.grid(idx_tol) = inf;
        times.scale(idx_tol) = inf;
        for i=1:rep
            pre_t = tic;
            pre = fse.precomp(opt);
            times.pre(idx_tol) = min(toc(pre_t), times.pre(idx_tol));
                        
            [uk, fourier_time_detail] = fse.fourier_sum(x, f, opt, pre);
            [ur, real_time_detail]    = fse.real_sum(x, f, opt);
            
            times.real(idx_tol) = min(real_time_detail.eval, ...
                                    times.real(idx_tol));
            times.fourier(idx_tol) = min(fourier_time_detail.total, ...
                                       times.fourier(idx_tol));
            times.fft(idx_tol) = min(fourier_time_detail.fft, ...
                                   times.fft(idx_tol));
            times.grid(idx_tol) = min(fourier_time_detail.grid+fourier_time_detail.int, ...
                                    times.grid(idx_tol));
            times.scale(idx_tol) = min(fourier_time_detail.scale, ...
                                     times.scale(idx_tol));            
        end
        us = fse.self(f, opt);
        ue = uk + ur + us;
        errors(idx_tol) = rms(ue - uref) / rms(uref)
    end
    times.fsewald = times.real + times.fourier + times.pre;

    plot(errors, times.fsewald, '.-', 'color',col{idx_kernel}, ...
         'DisplayName', [kernel ' FSE'])
    hold on
    errors
    grid on
    xlabel('rms error (rel.)')
    ylabel('time [s]')
end
legend('toggle')
set(gca,'xdir','reverse')
set(gca,'xscale','log')

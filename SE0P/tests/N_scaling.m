clear


kernels = {'stokeslet', 'stresslet', 'rotlet'};

Nlist = [600 1200 2500 5000:5000:40000];

rho = 2500;
Llist = (Nlist/rho).^(1/3);
col = {'b','r','k'};

rep = 6; % Repetitions

opt.xi = 7;
opt.oversampling = 1+sqrt(3);
opt.no_extra_support = true;

rms = @(u) sqrt(1/size(u,1)*sum(u(:).^2));


for idx_kernel = 1:numel(kernels)
    kernel = kernels{idx_kernel};
    fse = fsewald(kernel);
    
    if strcmp(kernel,'stokeslet')
        fse.real_sum = @stokeslet_real_space_rangesearch;
    end
    
           
    switch kernel
        % Optimized at N = 20000 for error < .51e-8
      case 'stokeslet'       
        MoverL = 48/2;
        opt.P = 16;
        opt.rc = 0.63;
        Mlist = [16  20  24  30  38  44  48  52 54 56 62]-2;
      case 'stresslet'
        MoverL = 50/2;
        opt.P = 16;
        opt.rc = 0.63;
        Mlist = [16  20  24  30  38  44  48  52 54 56 60];
      case 'rotlet'
        MoverL = 38/2;
        opt.P = 16;
        opt.rc = 0.58;
        Mlist = [12 16 20 24 30 36 38 42 42 46 48];
    end
    
    errors = [];
    times.direct = zeros(size(Nlist));
    for idx_N = 1:numel(Nlist);
        rng(1);
        L = Llist(idx_N);
        opt.M = round(MoverL*L*[1 1 1]);
        
        % just to modify some M values in stresslet and stokeslet
        opt.M = Mlist(idx_N)*[1 1 1];

        N = Nlist(idx_N);
        box = [L L L];
        opt.box = box;
        [x f] = fse.generator(N, box);
        
        times.direct(idx_N) = inf;
        times.pre(idx_N) = inf;
        times.real(idx_N) = inf;
        times.fourier(idx_N) = inf;
        times.fft(idx_N) = inf;
        times.grid(idx_N) = inf;
        times.scale(idx_N) = inf;
        for i=1:rep
            t = tic();
            uref = fse.direct_sum(x, f, box);
            times.direct(idx_N) = min(toc(t), times.direct(idx_N))

            pre_t = tic;
            pre = fse.precomp(opt);
            times.pre(idx_N) = min(toc(pre_t), times.pre(idx_N));
                        
            [uk, fourier_time_detail]= fse.fourier_sum(x, f, opt, pre);
            [ur, real_time_detail] = fse.real_sum(x, f, opt);
            
            times.real(idx_N) = min(real_time_detail.eval, ...
                                    times.real(idx_N));
            times.fourier(idx_N) = min(fourier_time_detail.total, ...
                                       times.fourier(idx_N));
            times.fft(idx_N) = min(fourier_time_detail.fft, ...
                                   times.fft(idx_N));
            times.grid(idx_N) = min(fourier_time_detail.grid+fourier_time_detail.int, ...
                                    times.grid(idx_N));
            times.scale(idx_N) = min(fourier_time_detail.scale, ...
                                     times.scale(idx_N));            
        end
        us = fse.self(f, opt);        
        ue = uk + ur + us;
        errors(idx_N) = rms(ue - uref) / rms(uref);
    end
    times.fsewald = times.real + times.fourier + times.pre;
    figure(1)
    plot(Nlist, times.direct, '--', 'color',col{idx_kernel}, ...
         'DisplayName', [kernel ' direct'])    
    hold on
    plot(Nlist, times.fsewald, '.-', 'color',col{idx_kernel}, ...
         'DisplayName', [kernel ' FSE'])
    errors
    grid on
    xlabel('N')
    ylabel('time [s]')
    
    if strcmp(kernel, 'stokeslet')
        figure(2)
        plot(Nlist, times.fsewald, 'k.-', 'DisplayName', ...
             [kernel ' FSE'])
        hold on
        plot(Nlist, times.fourier, 'b*-', 'DisplayName', ...
             [kernel ' Fourier sp.'])
        plot(Nlist, times.real, 'ro-', 'DisplayName', ...
             [kernel ' Real sp.'])
        plot(Nlist, times.pre, '^-', 'color',colorbox(4), 'DisplayName', ...
             [kernel ' Precomp.'])
        grid on
        xlabel('N')
        ylabel('time [s]')
        
        figure(3)
        plot(Nlist, times.fourier, 'k.-', 'DisplayName', ...
             [kernel ' Fourier sp.'])
        hold on
        plot(Nlist, times.grid, 'b*-', 'DisplayName', ...
             [kernel ' Grid'])
        plot(Nlist, times.fft, 'ro-', 'DisplayName', ...
             [kernel ' FFT'])
        plot(Nlist, times.scale, '^-', 'color',colorbox(4), 'DisplayName', ...
             [kernel ' Scale'])
        grid on
        xlabel('N')
        ylabel('time [s]')
    end

    times.fsewald_nopre = times.real + times.fourier;
    figure(4)
    plot(Nlist, times.direct, '--', 'color',col{idx_kernel}, ...
         'DisplayName', [kernel ' direct'])    
    hold on
    plot(Nlist, times.fsewald_nopre, '.-', 'color',col{idx_kernel}, ...
         'DisplayName', [kernel ' FSE'])
    errors
    grid on
    xlabel('N')
    ylabel('time [s]')
end
legend('toggle')

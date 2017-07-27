clear

M = 50;
L = 3;
N = 10000;

input('ENTER to overwrite current plots, Ctrl-c to cancel');

kernels = {'stokeslet','rotlet','stresslet'};
for i=1:3
    sfigure(i);
    d = show_fourier_error(kernels{i}, 'M',M, 'L',L, 'N',N);
    sfigure(i);    
    clf
    semilogy(d.Kxi, d.err_rms_rel, '.', 'DisplayName','Measured')    
    hold on
    semilogy(d.Kxi, d.est_rel, '-', 'DisplayName','Estimate')
    ylim([1e-15 1])
    xlim([0 15])
    grid on
    xlabel('$\kmax / \xi$','interpreter','latex')
    write_tikz(i, ['../latex/fig/fourier_err_' kernels{i}]);
end

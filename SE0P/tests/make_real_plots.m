clear

L = 3;
N = 2000;

input('ENTER to overwrite current plots, Ctrl-c to cancel');

kernels = {'stokeslet','rotlet','stresslet'};
for i=1:3
    sfigure(i);
    d = show_real_error(kernels{i}, 'L',L, 'N',N);    
    sfigure(i);    
    clf
    semilogy(d.xirc, d.err_rms_rel, '.', 'DisplayName','Measured')    
    hold on
    semilogy(d.xirc, d.est_rel, '-', 'DisplayName','Estimate')
    ylim([1e-16 1])
    xlim([0 7])
    grid on
    xlabel('$\xi r_c$','interpreter','latex')
    write_tikz(i, ['../latex/fig/real_err_' kernels{i}]);
end

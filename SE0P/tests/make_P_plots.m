
M = 40;
L = 1;
Kxi = 20;
sf = 4;
N = 1000;

input('ENTER to overwrite current plots, Ctrl-c to cancel');

for plotnum = [1 2]
    sfigure(plotnum);
    no_extra = (plotnum==2);
    d = show_P_error('stokeslet', 'M', M, 'N', N, 'Kxi', Kxi, 'oversampling', sf,...
                     'no_extra_support', no_extra);
    sfigure(plotnum);
    clf
    styles = {'o-k', '.-k', 's-k', '*-k', 's-k'};
    for j=1:numel(d.Plist)
        semilogy(d.Kxi, d.err_rms(j,:) / d.ref_max, styles{j}, ...
                 'DisplayName', sprintf('P=%d',d.Plist(j)))    
        hold on
    end
    ylim([1e-14 1])
    xlim([0 Kxi])
    grid on
    xlabel('$\kmax / \xi$','interpreter','latex')
    h = legend('toggle');
    set(h, 'Location','SouthWest')
    drawnow
    
    write_tikz(plotnum, ['../latex/fig/stokeslet_P_err_' num2str(plotnum)]);
end

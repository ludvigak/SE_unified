clear
fse_warnings('off');

L = 2;
M_min = 28;
M = 52;
N = 20000;
xi = 7;
Kxi = pi*M/xi/L;
%sf = 1 + sqrt(3);
sf = 3;
Plist = 4:2:22;

sfigure(1);
d1 = show_P_error('stokeslet', 'M', M, 'N', N, 'Kxi', Kxi, 'L', L, 'oversampling', sf,...
                 'no_extra_support', true, 'Plist', Plist, 'M_min', M_min);
d2 = show_P_error('stresslet', 'M', M, 'N', N, 'Kxi', Kxi, 'L', L, 'oversampling', sf,...
                 'no_extra_support', true, 'Plist', Plist, 'M_min', M_min);

markers = '.ox+*sdph<>^v';
sfigure(2);
clf
for j=1:numel(Plist)
    semilogy(d1.Mlist, d1.err_rms(j,:) / d1.full_ref_rms, [markers(j) '-r'], ...
             'DisplayName', sprintf('stokeslet, P=%d',Plist(j)))    
    hold on
end
for j=1:numel(Plist)
    semilogy(d2.Mlist, d2.err_rms(j,:) / d2.full_ref_rms, [markers(j) '-k'], ...
             'DisplayName', sprintf('stresslet, P=%d',Plist(j)))        
end
grid on
xlabel('M')

title(sprintf('N=%d, \\xi=%g, L=%d, s_f=%.2f',N,d1.xi,L,sf))

h = legend('toggle');
set(h, 'Location','EastOutside')
drawnow
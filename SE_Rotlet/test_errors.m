clear


rng(1);

s = 1;

L = 1*s;
M_max = 80/s;
xi = 28/s;

N = 500;
t = 1-2*rand(N,3);

box = L*[1 1 1];
x = L*rand(N,3);
xe = L*rand(N,3);


shells = 15;
k_shells = 15;


opt.P = 32;
opt.box = box;

legh = [];
leglist = {};
clf
hold all
co = get(gca,'ColorOrder');
set(gca,'ColorOrder',co(1:2,:))
set(gca,'LineStyleOrder','.-|o-|*-')
set(gca,'yscale','log');
for i=1:3
    opt.M = (M_max+20)*box;
    ref = SE_Rotlet(xe, x, t, xi, opt);
    ref_max = norm(ref(:), inf);

    Mlist = 2:2:M_max;
    %Mlist = 4:4:M_max;
    err_inf = [];
    err_rms = [];
    for i=1:numel(Mlist)
        this_opt = opt;
        this_opt.M = Mlist(i)*box;    
        uk = SE_Rotlet(xe, x, t, xi, this_opt);
        err = uk - ref;
        err_inf(i) = norm(err(:), inf) / ref_max;
        err_rms(i) = sqrt(1/N*sum(err(:).^2)) / ref_max;
    end

    K = pi*Mlist;
    legh(end+1) = semilogy(Mlist, err_rms)       
    
    Q = sum(t(:).^2);
    V = prod(box);

    % original force (Kolafa)
    %semilogy(K, sqrt(16*xi^2*Q./(pi*V*K)) .* exp(-K.^2/4/xi^2) / ref_max,'--')

    %semilogy(Mlist, sqrt(2/3)*sqrt(4*xi^2*Q./(pi*V*K)) .* exp(-K.^2/4/xi^2) / ref_max,'-')
    semilogy(Mlist, sqrt(8*xi^2*Q./(3*pi*V*K)) .* exp(-K.^2/4/xi^2) / ref_max,'-')
    
    ylim([1e-16 1])
    grid on
    xlabel('$M = K/\xi$','interpreter','latex')
    ylabel('$e_{rms}$ (rel.)','interpreter','latex')
    leglist{end+1} = sprintf('$\\xi=%g$',xi);
    legend(legh,leglist,'interpreter','latex');
    drawnow
    xi = xi/2;
end
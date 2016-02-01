clear

savedir = 'fig/';

cases = [1 2 3 4 5];
cases = 2;

if numel(cases)>1, close all, end

for testcase = cases

switch testcase
    case 1
        xilist = [10 20 30 40 50];
        N = 5000;
        Mlist = 3:2:125; 
        box = [1 1 1];
    case 2
        xilist = [5 10 15 20 25];
        N = 7000;
        Mlist = 3:2:135;
        box = [2 2 2];
        
    case 3
        xilist = [5 10 15 20 25 30 40];% 50];
        N = 7000;
        Mlist = 3:2:93;
        box = [1 2 2];
    case 4
        xilist = [8 10 15 20 30];
        N = 2000;
        Mlist = 3:2:81;
        box = [1 1 5];

    case 5
        xilist = [10 15 20];
        N = 1000;
        Mlist = 3:2:81;
        box = [2 4 6];
    otherwise
            error('Unknown test case')
end

generate_ref_state = false;
generate_ref_solutions = false;
generate_errors = false;

% ******************************

kmax = ( (Mlist-1)/2 );
boxstr = sprintf('box(%g,%g,%g)',box(1),box(2),box(3));
SE_opt.box = box;
ubox = box/min(box);

filename_state = ['tests/refmat/test_error_kmax_ref_state_' boxstr '.mat'];
if generate_ref_state
    disp('Generating reference state...')
    [x f nvec] = generate_state(N,box);
    idx = 1:N;
    save(filename_state,...
          'N','x','f','nvec','box','idx');
    disp('done.')
    fprintf('Saved %s\n',filename_state);
else
    load(filename_state);
end

if generate_ref_solutions
    for xi=xilist
        filename_sol = sprintf('tests/refmat/test_error_kmax_ref_sol_%s_xi=%.3f.mat',boxstr,xi);    
        SE_opt.P = 32;
        
        SE_opt.M = (2*Mlist(end))*ubox;
        SE_static  = SE_Stresslet_pre(x,xi,SE_opt);
        UREF2 = SE_Stresslet(idx,x,f,nvec,xi,SE_opt,SE_static);
        
        SE_opt.M = (2*Mlist(end)+4)*ubox;
        SE_static  = SE_Stresslet_pre(x,xi,SE_opt);
        UREF = SE_Stresslet(idx,x,f,nvec,xi,SE_opt,SE_static);
        
        disp('*************************')
        fprintf('************************* UREF precision: %g\n',norm(UREF(:)-UREF2(:),inf));
        disp('*************************')
        
        save(filename_sol,...
              'N','x','f','nvec','box','idx','UREF','SE_opt');
        fprintf('Saved %s\n',filename_sol);
    end
end



err = [];
for j=1:numel(xilist);
    xi = xilist(j);
    filename_sol = sprintf('tests/refmat/test_error_kmax_ref_sol_%s_xi=%.3f.mat',boxstr,xi);    
    load(filename_sol)    
    filename_errs = sprintf('tests/refmat/test_error_kmax_errs_%s_xi=%.3f.mat',boxstr,xi);  
    if generate_errors
        for i=1:numel(Mlist)
            if i>2
                err_update = abs(errxi(i-1)-errxi(i-2));
                if err_update < 1e-12
                    errxi(i) = errxi(i-1);
                    errxi_j(i,:) = errxi_j(i-1,:);
                    disp('Converged')
                    continue
                end
            end
            fprintf('M = %d, xi = %g\n', Mlist(i), xi);
            SE_opt.M = Mlist(i)*ubox;
            SE_opt.P = 32;
            SE_static  = SE_Stresslet_pre(x,xi,SE_opt);
            u = SE_Stresslet(idx,x,f,nvec,xi,SE_opt,SE_static);
            e = u - UREF;
            esq = e.^2; 
            erms = sqrt( sum(esq(:))/(3*N) );
            errxi(i) = erms;
            errxi_j(i,:) = sqrt( sum(esq,1)/(N) );
        end
        save(filename_errs,'errxi','errxi_j');
        fprintf('Saved %s\n',filename_errs);
    else
        load(filename_errs);
    end
    err = [err errxi(:)];
end

leglist = {};
for i=1:numel(xilist)
    leglist{i} = sprintf('\\xi=%g', xilist(i));
end

%% Plot
sfigure(1); clf, publication_fig
L = min(box);
K = 2*pi*kmax/L;
semilogy(K, err,'.-b')
titlestr=['Truncation error of k-space part (RMS), \xi = '...
        sprintf('%g,',xilist)];
title(titlestr(1:end-1));
xlabel('K')
ylabel('E_{RMS}')
grid on
hold on
set(gca,'Box','on')
% legend(leglist,'Location','BestOutside')

% ESTIMATE
S2 = zeros(3,3);
for l=1:3
    for m=1:3
        S2(l,m) = sum( ( f(:,l).*nvec(:,m) ).^2 );
    end
end
Q = sum(S2(:));
V = prod(box);
L = min(box);

% Create estimate
est = zeros(numel(xilist),numel(kmax));
for j=1:numel(xilist);
    xi = xilist(j);        
    kc = 2*kmax*pi/L;
    kcxi = kc/xi;

    % Simple
    est_simple = sqrt(Q/V)*sqrt(L)*xi^1*exp(-kcxi.^2/4).*exp(.43*kcxi);
    
    est(j,:) = est_simple;
end

semilogy(K, est,'--r')
axis tight
ylim([1e-12 200])
publication_fig
% Figure 1 done

% Plot accuracy of estimate
sfigure(2); clf,publication_fig
plot(kc'*xilist.^-1, est'./err, '.-', [0 15],[1 1],'--k')
ylabel('est / err')
xlabel('K / \xi')
ylabel('Estimate / Measured (RMS)')
title('Accuracy of k-space truncation error estimate')
xlim([4 10])
grid on
ylim([0.8 2])
legend(leglist,'Location','NorthWest')
publication_fig

% Figure 2 done

if (0)
    write_fig(1, [savedir 'e_rms_k'])
    write_fig(2, [savedir 'e_rms_k_acc'])
end


% % Attempt data collapse
% sfigure(3);clf
% x = 2*pi/L*kmax'*(xilist).^-1;
% y = err;
% y(err<1e-10) = 0;
% y = y.*exp((x/2).^2);
% % y = y.*exp(-0.43*x);
% y = y.*exp(-pi/7*x);
% % y = y.*x.^(1/7);
% y = bsxfun(@times,y,xilist.^-1); % 1 is more right!
% semilogy(x,y,'.-')
% xlabel('k_c / \xi')
% 
% xmin = min(x(y>0));
% xmax = max(x(y>0));
% xlim([xmin xmax])
% ylim([min(y(:)) max(y(:))])
% 
% y(y==0) = NaN;
% xinterp = xmin:0.001:xmax;
% yinterp = [];
% 
% for i=1:numel(xilist)
%     yinterp = [yinterp; interp1(x(:,i), y(:,i), xinterp,'pchip', inf)];
% end
% hold on
% grid on
% 
% % Relation between collapsed curves
% sfigure(4); clf
% xinterp = xinterp';
% yinterp = yinterp';
% for i=1:size(yinterp,2)
%     plot(xinterp,yinterp(:,i)./yinterp(:,end-2));
%     hold all
% end
% xlabel('k_c / \xi')
% ylim([0 2])
% grid on



% Complete collapse in figure 5
sfigure(5)
x = 2*pi/L*kmax'*(xilist).^-1;
y = err;
y = bsxfun(@times,y,xilist.^-1);
y = y/sqrt(Q/V*L);
semilogy(x,y,'.k'), hold on
xlim([0 15])


end
%%

if numel(cases)>1
    % Do something with data from complete collapse
    figure(5)
    lh = findall(gcf,'type','line');

    xl = [];
    yl = [];
    for i=1:numel(lh)
        h=lh(i);
        xl = [xl get(h,'xdata')];
        yl = [yl get(h,'ydata')];
    end
    i = xl<15;
    xl = xl(i);
    yl = yl(i);

    [xl,p] = sort(xl);
    yl = yl(p);
    disp('Saving data collapse')
    save('tests/refmat/kmax_err_allcollapse.mat','xl','yl')
end
    
return
%%
load('tests/refmat/kmax_err_allcollapse.mat')
figure(7), clf, publication_fig

ifit = 0<=xl & xl<=110;
xfit = xl(ifit);
yfit = yl(ifit);
% yfit = log(yfit)+(xfit/2).^2;

k = .43;
% k = exp(1)/2/pi;
m = 2.59-6*k;
m=0;
% plot(xfit,yfit,'.k', [0 xfit],k*[0 xfit]+m,'r')

semilogy(xfit,yfit.*exp((xfit/2).^2),'.k', [0 xfit],exp(k*[0 xfit]),'r')
xlabel('x')
xlim([4 12])
grid on
legend('g(x) \cdot e^{x^2/4}',...
        sprintf('e^{%.2fx}',k),'Location','Best')
publication_fig

figure(6), clf, publication_fig
semilogy(xl, yl, '.k', x, exp(-(x/2).^2+k*x)*exp(m),'r')
xlabel('x')
xlim([4 15])
ylim([1e-15 1])
legend('g(x)',...
        sprintf('e^{-x^2/4+%.2fx}',k),'Location','Best')
grid on

publication_fig

if (0)
    write_fig(6, [savedir 'e_rms_k_coll1'])
    write_fig(7, [savedir 'e_rms_k_coll2'])
end
return
% Spectral Ewald, basic accuracy/convergence computation
% Dag Lindbo, dag@kth.se

clear all,  close all

box = [1 1 1]; % domain
N = 10;        % numbver of charged particles
xi = 4;        % ewald parameter

P = 5:2:25; % support 
m = [5 6.5 8]; % shape

% grid
SE_opt.M = 28*box
SE_opt.box = box;

% charge-neutral system
[x q] = SE_charged_system(N,box,'scalar');

% parameters for (reference) direct Ewald sum
ED_opt.layers = (SE_opt.M(3)-1)/2;
ED_opt.xi = xi;
ED_opt.box = box;

% compute FD Ewald sum
ref = SE3P_direct_fd_mex(1:N,x,q,ED_opt);

for i = 1:length(m)
    % for Gaussian shaper parameter, m
    
    SE_opt.m = m(i);
    leg{i} = sprintf('m = %.1f',SE_opt.m) % plot legend string (below)
    for j = 1:length(P)
        % for Gaussian support, P
        
        SE_opt.P = P(j);
        u = spectral_ewald(1:N,x,q,xi,SE_opt);

        % compute RMS error
        e = (u - ref).^2; 
        err(i,j) = sqrt( sum(e)/N );
        
    end
end

% same thing for m = m(P)
SE_opt = rmfield(SE_opt,'m')
leg{i+1} = sprintf('m(P)')
for j = 1:length(P)
    SE_opt.P = P(j);
    u = spectral_ewald(1:N,x,q,xi,SE_opt);
    e = (u - ref).^2; 
    err(i+1,j) = sqrt( sum(e)/N );
end

sty = {'.-','*-','+-','r'}
hold on
for i = 1:size(err,1)
    plot(P,err(i,:),sty{i})
end

publication_fig
set(gca,'YScale','log')
set(gca,'YTick',[1e-15 1e-10 1e-5 1e-0])
xlabel('P')
ylabel('e_{rms}')
grid on
axis([P(1) P(end) 1e-16 1e2])
legend(leg,'Location','Best')
fname = sprintf('output/SE_accuracy_xi%d_M%d',xi,SE_opt.M(1));
write_fig(1,fname);

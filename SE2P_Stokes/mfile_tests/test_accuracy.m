% 2P Stokeslet Spectral Ewald, basic accuracy/convergence computation
% Dag Lindbo, dag@kth.se

clear all,  close all

SE_opt.box = [1 1 1]; % domain
N = 10;               % numbver of charged particles
xi = 4;               % ewald parameter

P = 5:2:25; % support 
s = 1:3;

% grid
SE_opt.M = 31;

% charge-neutral system
[x f] = SE_charged_system(N,SE_opt.box,'vector');

% parameters for (reference) direct Ewald sum
ED_opt.layers = (SE_opt.M-1)/2;
ED_opt.xi = xi;
ED_opt.box = SE_opt.box;

% compute FD Ewald sum
ref = SE2P_Stokes_direct_fd_mex(1:N,x,f,ED_opt);

for i = 1:length(s)
    % for oversampling factor, s
    SE_opt.sampling_factor=s(i);
    leg{i} = sprintf('s_f = %d',2*s(i)) % plot legend string (below)
    for j = 1:length(P)
        % for Gaussian support, P
        
        SE_opt.P = P(j);
        u = SE2P_Stokes(1:N,x,f,xi,SE_opt);

        % compute RMS error (first component)
        e = (u - ref).^2; 
        err(i,j) = sqrt( sum(e(:,1))/N );
        
    end
end

sty = {'.-','*-','+-'};
hold on
for i = 1:size(err,1)
    plot(P,err(i,:),sty{i})
end

set(gca,'YScale','log')
set(gca,'YTick',[1e-15 1e-10 1e-5 1e-0])
xlabel('P')
ylabel('e_{rms}')
grid on
axis([P(1) P(end) 1e-16 1e2])
legend(leg,'Location','Best')

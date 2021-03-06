
% Stresslet Spectral Ewald, basic accuracy/convergence computation

clear all

box = [1 1 1]; % domain
N = 10;        % numbver of charged particles
xi = 4;        % ewald parameter

eval_idx = 1:N; % eval points

P = 5:2:25; % support 
m = [5 6.5 8]; % shape

% grid
SE_opt.M = 31*box
SE_opt.box = box;

% charge-neutral system
[x f nvec] = generate_state(N,box);

% parameters for (reference) direct Ewald sum
ED_opt.layers = (SE_opt.M(1)-1)/2;
ED_opt.xi = xi;
ED_opt.box = box;

% compute FD Ewald sum
ref = stresslet_direct_fd( eval_idx, x, f, nvec, ...
        ED_opt.xi, ED_opt.box, ED_opt.layers);

for i = 1:length(m)
    % for Gaussian shaper parameter, m
    
    SE_opt.m = m(i);
    leg{i} = sprintf('m = %.1f',SE_opt.m) % plot legend string (below)
    for j = 1:length(P)
        % for Gaussian support, P
        
        SE_opt.P = P(j);
        u = SE_Stresslet(eval_idx,x,f, nvec,xi,SE_opt);

        % compute RMS error (first vector component)
        e = (u - ref).^2; 
        err(i,j) = sqrt( sum(e(:,1))/N );
        
    end
end

% same thing for m = m(P)
SE_opt = rmfield(SE_opt,'m')
leg{i+1} = sprintf('m(P)')
for j = 1:length(P)
    SE_opt.P = P(j);
    u = SE_Stresslet(eval_idx,x,f, nvec,xi,SE_opt);
    e = (u - ref).^2; 
    err(i+1,j) = sqrt( sum(e(:,1))/N );
end


%%
figure(3);
sty = {'.-','*-','+-','r'}
hold on
for i = 1:size(err,1)
    plot(P,err(i,:),sty{i})
end

% publication_fig
set(gca,'YScale','log')
set(gca,'YTick',[1e-15 1e-10 1e-5 1e-0])
xlabel('P')
ylabel('e_{rms}')
grid on
axis([P(1) P(end) 1e-16 1e2])
legend(leg,'Location','Best')
fname = sprintf('output/SE_accuracy_xi%d_M%d',xi,SE_opt.M(1));

status = 1;
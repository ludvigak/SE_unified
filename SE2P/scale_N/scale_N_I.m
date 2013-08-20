clear all, close all

epsi=1e-10;
SE_opt.P = 16;
SE_opt.sampling_factor = 4;

%epsi=1e-6;
%SE_opt.P = 10;
%SE_opt.sampling_factor = 3;

NN = round(linspace(1e3,1e5,10))
[xi M rc] = scale_params_rms(NN,epsi);
SE_opt.box = [1 1 1];
M = max(M,SE_opt.P);

% nicer transform lengths
M = [24    56    72    80    88    96   100   104   112   120]
%M = [14    34    44    48    56    60    64    68    72    72]

for j = 1:length(NN)
   
    N = NN(j);

    SE_opt.M = M(j);

    h = 1/M(j);
    P = SE_opt.P;
    SE_opt.zLim = [-h*P/2 1+h*P/2];
    
    [x q] = SE_charged_system(N,SE_opt.box,'scalar');
    SE_fgg_static = SE_FGG_precomp(x,xi(j),SE_opt);
    tic
    phi = spectral_ewald_2P(1:N,x,q,xi(j),SE_opt,SE_fgg_static);
    wtime(j) = toc;
end

c = [ones(size(NN')) NN']\wtime';
Nc = [NN(2) NN(end-1)]
plot(NN,wtime,'*', Nc, c(2)*Nc-0.05,'r')
%axis([0 NN(end) 0 1])

publication_fig
xlabel('N'), ylabel('time (s)'), 
legend('measured','linear','Location','NorthWest')
fname = sprintf('scale_N_rms_wtime_N%d_P%d',N(end),P)
write_fig(1,fname)

dat = [NN' rc' xi' M'];
fprintf('%d & %.2f & %.2f & %d \\\\ \n',dat')
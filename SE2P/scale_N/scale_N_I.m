clear all, close all

epsi=1e-10;
opt.P = 16;
opt.s = 4;
opt.s0= 2;

%epsi=1e-6;
%opt.P = 10;
%opt.s = 3;
%opt.s0= 2;

NN = round(linspace(1e3,1e5,10))
[xi M rc] = scale_params_rms(NN,epsi);
opt.box = [1 1 1];
M = max(M,opt.P);

% nicer transform lengths
M = [24    56    72    80    88    96   100   104   112   120]
%M = [14    34    44    48    56    60    64    68    72    72]

for j = 1:length(NN)
   
    N = NN(j);

    opt.M = M(j);
    opt.xi = xi(j);

    h = 1/M(j);
    P = opt.P;
    
    [x q] = vector_system(N,opt.box);
    tic
    phi = se2p_fourier_space(x,q,opt);
    wtime(j) = toc;
end

c = [ones(size(NN')) NN']\wtime';
Nc = [NN(2) NN(end-1)]
plot(NN,wtime,'*', Nc, c(2)*Nc,'r')
%axis([0 NN(end) 0 1])

publication_fig
xlabel('N'), ylabel('time (s)'), 
legend('measured','linear','Location','NorthWest')
fname = sprintf('scale_N_rms_wtime_N%d_P%d',N(end),P)
%write_fig(1,fname)

dat = [NN' rc' xi' M'];
fprintf('%d & %.2f & %.2f & %d \\\\ \n',dat')

% 1P Spectral Ewald + k0 term, basic accuracy/convergence computation
% Periodic in z and free in x and y.
% Davoud Saffar Shamshirgar, davoudss@kth.se
% 2016/04/12

clear all;
tol = 1e-6;

NN = [300 300 2400 19200 153600 1228800];
LL = [10 10 20 40 80 160];
MM = [16 16 32 64 128 256];

% the followings are valid for wall systems
NL = [2 2 3 5 9 11];
SL = [2 2 2 2 2 2];
S0 = [2 2 2.5 2.5 2.5 2.5];
%rel_error_1p= [2.81e-6 2.81e-6 1.26e-5 2.84e-5 1.57e-5 8.93e-6];

opt.xi = .7186;        % ewald parameter
opt.P = 12;

reps = 5;          % number of repetitions
compute_ref = 1;
compute_error = 1;

for i=1:4%numel(NN)
    N = NN(i);
    L = LL(i);  
    M = MM(i);
    opt.box = [L L L];
    
    fp = fopen(sprintf('wall%d.bin',N),'r');
    % fp = fopen('wall300.bin','r');
     x=fread(fp,'double');
     fclose(fp);
     x = reshape(x,N,3);
     q = (-1).^(1:N)';
    
    % compute FD Ewald sum
    if(N<=10 && compute_ref)
        % parameters for (reference) direct Ewald sum
        opt.layers = (M/2+4)*[1 1 1];
        ref1p = se1p_k0_direct_force(1:N,x,q,opt,true);
        ref1p = ref1p + se1p_fourier_space_direct_force(1:N,x,q,opt,true);

        ref3p = SE3P_direct_fd_force_mex(1:N,x,q,opt);
    end
    
    fp = fopen(sprintf('ref1p_%d_force.bin',N),'r');
    ref1p = fread(fp,[N 3],'double');
    fclose(fp);    
    
    % for oversampling factor, s over the whole grid
    opt.s0 = S0(i);
    opt.sl = SL(i);
    opt.nl = NL(i);   
    opt.M = M;
   
    [times.pre, times.grid, times.fft, times.scale, times.int] = deal(inf);
    
    for k=1:reps
        [u1p,wt] = se1p_fourier_space_force(x, q, opt);
        times = time_min(times,wt);        
    end
    times_1p(i) = times;
    total_1p(i) = sum(struct2array(times));
    
    opt.M = M*[1 1 1];
    [times.pre, times.grid, times.fft, times.scale, times.int] = deal(inf);
    for k = 1:reps
        [u3p,wt] = se3p_fourier_space_force(x, q, opt);
        times = time_min(times,wt);
    end
    times_3p(i) = times;
    total_3p(i) = sum(struct2array(times));
    
    if(~compute_ref && compute_error)
        rms_1p = rms(ref1p-u1p)/rms(ref1p)
        rms_3p = rms(ref3p-u3p)/rms(ref3p)
    end
end

total_1p(1) = [];
total_3p(1) = [];

figure(1)
hold on

% SE
NN = [300 2400 19200 153600 1228800];
NN = NN(1:numel(total_1p));
plot(NN,total_1p,'s-','markersize',4,'color',colorbox(4),'DisplayName','SE1P')
plot(NN,total_3p,'ro-','markersize',4,'DisplayName','SE3P')

% NUFFT
NN = [300 2400 19200 153600 1228800];
nestler1p=[4e-2 2.3e-1 2.2e0 1.3e1 6.7e2]';
nestler3p = [1e-2 7e-2 5.4e-1 4.48e0 3.6e1]';
plot(NN,nestler1p','k^-','markersize',4,'DisplayName','NFFT1P')
plot(NN,nestler3p','k*-','markersize',4,'DisplayName','NFFT3P')

legend toggle
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'xtick',10.^(2:6))
set(gca,'ytick',10.^(-3:3))
axis([1e2 3e6 1e-4 1e3])
grid on
box on
ylabel('time [s]')
xlabel('N')


figure(2)
hold on
% SE
NN = [300 2400 19200 153600 1228800];
NN = NN(1:numel(total_1p));
plot(NN,total_1p./total_3p,'bs-','markersize',4,'DisplayName','SE1P/SE3P')
hold on

% NUFFT
NN = [300 2400 19200 153600 1228800];
plot(NN,nestler1p./nestler3p,'r^-','markersize',4,'DisplayName','NFFT1P/NFFT3P')

legend toggle
set(gca,'xscale','log')
%set(gca,'yscale','log')
set(gca,'xtick',10.^(2:6))
%set(gca,'ytick',10.^(-3:3))
%axis([1e2 3e6 1e-3 1e3])
grid on
box on
ylabel('time [s]')
xlabel('N')

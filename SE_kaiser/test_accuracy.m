clf

rng(1)
box = [1 1 1];   % domain
N = 1000;         % number of charged particles

M0 = 28;

opt.M = M0*box;
opt.xi = pi*M0 / 12;opt.xi=4;
opt.rc = 6 / opt.xi;
opt.box = box;

opt.layers = (opt.M(3)-1)/2;

opt.beta = 2*pi*.98^2;
p = [0.0059524 -0.19643 6.1905];
opt.beta = polyval(p,opt.xi);

% charge-neutral system
[x q] = SE_charged_system(N,box,'scalar');

idx = 1:10;
% compute FD Ewald sum
ref = SE3P_direct_fd_mex(idx,x,q,opt);

rms = @(x) sqrt(sum(x.^2)/numel(x));

Plist = 4:2:28;
[EK EG TKgrid TGgrid TKfft TGfft TK TG] = deal([]);
for P=Plist
    opt.P = P;
    [u time u0 time0]= spectral_ewald(1:N,x,q,opt);
    EK(end+1)    = rms(u(idx)-ref)/rms(ref);
    EG(end+1)    = rms(u0(idx)-ref)/rms(ref);
    TKgrid(end+1)= time.grid+time.int;
    TGgrid(end+1)= time0.grid+time0.int;
    TKfft(end+1) = time.fft;
    TGfft(end+1) = time0.fft;
    TK(end+1)    = sum(time.grid+time.int+time.fft+time.scale);
    TG(end+1)    = sum(time0.grid+time0.int+time0.fft+time0.scale);
end

figure(1)
semilogy(Plist,EK,'b.-',Plist,EG,'ro-')
hold on
semilogy(Plist,Plist.^(-5).*exp(-opt.beta/50*Plist.^2), ...
	'b:',Plist,exp(-.98^2*pi*Plist/2),'r--')
legend('Kaiser','Gaussian')
xlabel('P')
ylabel('e_{rms} rel.')
ylim([1e-15 1e0])
return
figure(2)
plot(Plist,TKgrid,'b.-',Plist,TGgrid,'ro-')
hold on
plot(Plist,TKfft,'b--',Plist,TGfft,'r--')
xlabel('P')
ylabel('time [s]')

figure(3)
semilogx(EK,TK,'b.-',EG,TG,'ro-')
legend('Kaiser','Gaussian')
ylabel('time [s]')
xlabel('e_{rms} rel.')


%fprintf('Kaiser  : %g\n',rms(u(idx)-ref)/rms(ref))
%fprintf('Gaussian: %g\n',rms(u0(idx)-ref)/rms(ref))

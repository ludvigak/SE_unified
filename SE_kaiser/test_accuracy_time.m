%clear all

rng(1)
box = [1 1 1];   % domain
N = 2000000;         % number of charged particles
M0 = 36;

opt.M = M0*box;
opt.xi = pi*M0 / 12;opt.xi=4;
opt.rc = 6 / opt.xi;
opt.box = box;

opt.layers = (opt.M(3)-1)/2;

opt.beta = 2*pi*.98^2;
p = [0.0059524 -0.19643 6.1905];
opt.beta = polyval(p,opt.xi);opt.beta = 2*pi*.97^2;

% charge-neutral system
[x q] = SE_charged_system(N,box,'scalar');
opt.P = 20;

idx = 1:N;
%ref= spectral_ewald_gaussian(idx,x,q,opt);
ref= spectral_ewald_kaiser(1:N,x,q,opt);

% compute FD Ewald sum
%idx = 1:10;
%ref = SE3P_direct_fd_mex(idx,x,q,opt);

opt.P = 16;
[u0, time0]= spectral_ewald_gaussian(1:N,x,q,opt);

opt.P = 10;
[u, time]= spectral_ewald_kaiser(1:N,x,q,opt);

EK    = rms(u(idx)-ref)/rms(ref);
EG    = rms(u0(idx)-ref)/rms(ref);

TKgrid= time.grid;
TKint = time.int;
TKfft = time.fft+time.ifft;
TKpre = sum(time.grid+time.int+time.fft+time.ifft+time.scale+time.pre+time.prefft);
TK    = sum(time.grid+time.int+time.fft+time.ifft+time.scale+time.pre);

TGgrid= time0.grid;
TGint = time0.int;
TGfft = time0.fft+time0.ifft;
TGpre=sum(time0.grid+time0.int+time0.fft+time0.ifft+time0.scale+time0.pre+time0.pregauss);
TG  =sum(time0.grid+time0.int+time0.fft+time0.ifft+time0.scale+time0.pre);

var={'Function','Error','grid','int','fft','TOTpre','TOTnopre'};
table({'KAISER','GAUSSIAN'}',[EK EG]',[TKgrid TGgrid]',[TKint TGint]',[TKfft TGfft]',[TKpre TGpre]',[TK TG]','variablenames',var)



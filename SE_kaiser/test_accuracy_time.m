clear all
rng(1)
box = [1 1 1];   % domain
N = 500000;         % number of charged particles

M0 = 28;

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
opt.P = 16;

idx = 1:N;
ref= spectral_ewald_gaussian(idx,x,q,opt);

% compute FD Ewald sum
%idx = 1:10;
%ref = SE3P_direct_fd_mex(idx,x,q,opt);

opt.P = 8;
[u0, time0]= spectral_ewald_gaussian(1:N,x,q,opt);

opt.P = 6;
[u, time]= spectral_ewald_kaiser(1:N,x,q,opt);

EK    = rms(u(idx)-ref)/rms(ref);
EG    = rms(u0(idx)-ref)/rms(ref);

TKgrid= time.grid;
TKint = time.int;
TKfft = time.fft+time.ifft;
TK    = sum(time.grid+time.int+time.fft+time.ifft+time.scale+time.pre+time.prefft);

TGgrid= time0.grid;
TGint = time0.int;
TGfft = time0.fft+time0.ifft;
TG    = sum(time0.grid+time0.int+time0.fft+time0.ifft+time0.scale+time0.pre);

var={'Function','Error','grid','int','fft','total'};
table(['K' 'G']',[EK EG]',[TKgrid TGgrid]',[TKint TGint]',[TKfft TGfft]',[TK TG]','variablenames',var)



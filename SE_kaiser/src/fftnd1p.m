function [x,xr,x0] = fftnd1p(x,opt)
%FFTND 3-dimensional discrete Fourier Transform.
%   FFTND(X,OPT) returns 3-dimensional DFT of 
%   the 3D array X with size M with two periodic directions
%   and 1 free directions in X vector and XR and X0. XR contains the LOCAL_PAD
%   result and the corresponding values in X vector are not updated. S is 
%   the oversampling factor on LOCAL_PAD and S0 is the oversampling factor 
%   on zero mod, all in free directions.
    
% Note: For smplicity, LOCAL_PAD contains zero mode as well but it
% will be updated by the zero mode later on.

%
%   INPUT:
%       X:          Input vector
%       Mx,My,Mz:   Size of the input vector (Mx,My,Mz)
%       LOCAL_PAD:  List of modes to apply a large oversampling (LOCAL_PAD=3:5).
%       S0:         Oversampling factor on zero mode
%       S:          Oversampling factor on LOCAL_PAD.
%                   (usally 2).
%
%   OUTPUT:
%       X:          Output vector of size (M,M,M)
%       XR:         Output vector of size (LOCAL_PAD,LOCAL_PAD,SR*M);
%       X0:         Output vector of size (1,1,S0*M);


Fx = fft(x ,opt.M);	       % 1D fft in x
x  = fft(fft(Fx,opt.My,2),opt.Mz,3);    % 2D fft with no padding in y,z
xr = fft(fft(Fx(opt.local_pad,:,:),round(opt.My*opt.s),2),round(opt.Mz*opt.s),3);% 1D fft with padding in z

Fx=squeeze(Fx(1,:,:));
x0 = fft2(Fx,round(opt.My*opt.s0),round(opt.Mz*opt.s0));

end

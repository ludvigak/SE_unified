function [x,xr,x0] = fftnd2p(x,opt)
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


Fxy = fft2(x ,opt.M, opt.M);	  % 2D fft in x and y
x   = fft(Fxy,opt.Mz,3);                   % 1D fft with no padding in z
xr  = fft(Fxy(opt.local_pad,opt.local_pad,:),round(opt.Mz*opt.s),3);% 1D fft with padding in z

F=Fxy(1,1,:);F=F(:);
x0 = fft(F,round(opt.Mz*opt.s0));
%x0 = transpose(x0);

end

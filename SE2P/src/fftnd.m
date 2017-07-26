function [x,x0,varargout] = fftnd(x,Mx,My,Mz,s0,s)
%FFTND 3-dimensional discrete Fourier Transform.
%   FFTND(X,M,S0,S,LOCAL_PAD,PER) returns 3-dimensional DFT of 
%   the 3D array X with size M with one periodic direction PER 
%   and 2 free directions in X vector and XRES. XRES contains the LOCAL_PAD
%   result and the corresponding values in X vector are not updated. S is 
%   the oversampling factor on LOCAL_PAD and S0 is the oversampling factor 
%   on the rest, all in free directions.
%
%   INPUT:
%       X:          Input vector
%       Mx,My,Mz:   Size of the input vector (Mx,My,Mz)
%       LOCAL_PAD:  List of modeshe to apply a large oversampling (LOCAL_PAD=3:5).
%       S0:         Oversampling factor on LOCAL_PAD
%       S:          Oversampling factor on the rest of the free domain.
%                   (usally 2).
%       PER:        Periodic diretion
%
%   OUTPUT:
%       X:          Output vector of size (S*M,S*M,M)
%       XRES:       Output vector of size (S0*M,S0*M,LOCAL_PAD);


Fxy = fft2(x ,Mx, My);                              % 2D fft in x and y
x   = fft(Fxy,round(Mz*s),3);                          % 1D fft with padding in z

x0 = fft(Fxy(1,1,:),Mz*s0);

end
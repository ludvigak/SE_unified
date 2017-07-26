function [x, varargout] = ifftnd(x,x0,Mx,My,Mz,s0,s)
%IFFTND 3-dimensional inverse discrete Fourier Transform.
%   IFFTND(X,M,S0,S,LOCAL_PAD,PER) returns 3-dimensional inverse DFT of 
%   the 3D array X with size M with one periodic direction PER 
%   and 2 free directions. S is the oversampling factor on LOCAL_PAD
%   and S0 is the oversampling factor on the rest, all in free directions.
%
%   INPUT:
%       X:          Input vector of size (S*M,S*M,M).
%       XRES:       LOCAL_PAD vector with S0 oversampling factor with size 
%                       (S0*M,S0*M,LOCAL_PAD).
%       M:          Size of the input vector in the periodic direction.
%       LOCAL_PAD:  List of modes to apply a large oversampling (LOCAL_PAD=3:5).
%       S0:         Oversampling factor on LOCAL_PAD.
%       S:          Oversampling factor on the rest of the free domain.
%                   (usally 2).
%       PER:        Periodic diretion.
%
%   OUTPUT:
%       X:          Output vector of size (M,M,M)


F0= ifft(x0,Mz*s0);

F1 = ifft(x,round(Mz*s),3);         % 1D fft on z

Fz = F1(:,:,1:Mz);		% update old values and restrict
Fz(1,1,:) = F0(1:Mz);           % restrict to (1,1,Mz)

x = ifft2(Fz,Mx,My);                              % 1D fft in z


end
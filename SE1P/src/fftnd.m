function [x,xres,x0mod,varargout] = fftnd(x,Mx,My,Mz,sl,s0,s,local_pad, zeromod, per)
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


% periodic direction should be z!
if per==1
    x = permute(x,[3 1 2]);
elseif per==2
    x = permute(x,[1 3 2]);
end
time = 0;

ft_t = tic;
Fz = fft(x  ,Mz ,3);                              % 1D fft in z
x  = fft2(Fz,round(Mx*s), round(My*s));           % This is also possible

if(~isempty(local_pad))
  xres = fft2(Fz(:,:,local_pad),Mx*sl,My*sl);% 2D fft local_pad with S0 oversampling 
else
    xres = 0;
end

x0mod = fft2(Fz(:,:,zeromod),Mx*s0,My*s0);          % 2D fft local_pad with S0 oversampling. 4 times oversampling is always needed.
time = toc(ft_t);
% switch back the dimensions
if per==1
    x = permute(x,[3 1 2]);
    xres = permute(xres,[3 1 2]);
    x0mod = permute(x0mod,[3 1 2]);
elseif per==2
    x = permute(x,[1 3 2]);
    xres = permute(xres,[1 3 2]);
    x0mod = permute(x0mod,[1 3 2]);
end

if(nargout==4)
    varargout{1} = time;
end
end
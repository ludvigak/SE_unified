function x = ifftnd2p(x,xr,x0,opt)
%IFFTND 3-dimensional inverse discrete Fourier Transform.
%   IFFTND(X,OPT) returns 3-dimensional inverse DFT of 
%   the 3D array X with size M with two periodic directions
%   and 1 free directions. S is the oversampling factor on LOCAL_PAD
%   and S0 is the oversampling factor on the zero mod, all in free directions.
%
% Note: For smplicity, LOCAL_PAD contains zero mode as well but it
% will be updated by the zero mode later on.

%
%   INPUT:
%       X:          Input vector of size (M,M,M).
%       XR:         LOCAL_PAD vector with S0 oversampling factor with size 
%                       (LOCAL_PAD,LOCAL_PAD,SR*M).
%       M:          Size of the input vector in the periodic direction.
%       LOCAL_PAD:  List of modes to apply a large oversampling (LOCAL_PAD=3:5).
%       S:          Oversampling factor on LOCAL_PAD.
%       S0:         Oversampling factor on the zero mod.
%
%   OUTPUT:
%       X:          Output vector of size (M,M,M)

F0= ifft(x0,round(opt.Mz*opt.s0)); % since this is a vector we skip 3

Fr= ifft(xr,round(opt.Mz*opt.s),3);

F1 = ifft(x,opt.Mz,3);		% 1D ifft on z

Fz = zeros(opt.M,opt.M,opt.Mz);

Fz = F1(:,:,1:opt.Mz);		% update old values and restrict
Fz(opt.local_pad,opt.local_pad,:) = Fr(:,:,1:opt.Mz);

Fz(1,1,:) = F0(1:opt.Mz);       % restrict to (1,1,Mz)

x = ifft2(Fz,opt.M,opt.M);      % 2D ifft in x and y

end

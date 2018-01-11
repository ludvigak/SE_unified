function x = ifftnd1p(x,xr,x0,opt)
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

global_pad= setdiff(1:opt.M,[opt.local_pad 1]);                 % global pad vector

F0= ifft2(x0,round(opt.My*opt.s0),round(opt.Mz*opt.s0)); 

Fr= ifft(ifft(xr,round(opt.My*opt.s),2),round(opt.Mz*opt.s),3);

F = ifft(ifft(x,opt.My,2),opt.Mz,3);		% 2D ifft on y and z

Fxy(opt.local_pad,:,:) = Fr(:,1:opt.My,1:opt.Mz);
Fxy(1,:,:)             = F0(1:opt.My,1:opt.Mz);
Fxy(global_pad,:,:)    = F(global_pad,1:opt.My,1:opt.Mz);

x = ifft(Fxy,opt.M);      % 1D ifft in x

end

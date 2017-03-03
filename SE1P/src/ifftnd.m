function [x, varargout] = ifftnd(x,xres,x0mod,Mx,My,Mz,sl,s0,s,local_pad, zeromod, per)
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


global_pad= setdiff(1:Mz,[local_pad zeromod]);                 % global pad vector

% periodic direction should be z!
if per==1
    x = permute(x,[3 1 2]);
    xres = permute(xres,[3 1 2]);
    x0mod = permute(x0mod,[3 1 2]);
elseif per==2
    x = permute(x,[1 3 2]);
    xres = permute(xres,[1 3 2]);
    x0mod = permute(x0mod,[1 3 2]);
end
time = 0;
ift_t = tic;
Fzres = ifft2(xres,Mx*sl,My*sl);               % 2D fft on local_pad in x and y
Fz0mod= ifft2(x0mod,Mx*s0,My*s0);              % 2D fft on local_pad in x and y, This is always 4 times oversampling
time = time + toc(ift_t);
if(~isempty(local_pad))
    Fz(:,:,local_pad) = Fzres(1:Mx,1:My,:);    % restrict to (Mx,My,local_pad)
end
Fz(:,:,zeromod) = Fz0mod(1:Mx,1:My);           % restrict to (Mx,My,zeromod)
ift_t = tic;
F1 = ifft2(x,round(Mx*s),round(My*s));         % 2D fft on the rest on x and y
time = time + toc(ift_t);
Fz(:,:,global_pad) = F1(1:Mx,1:My,global_pad);  % update old values and restrict
ift_t = tic;
x = ifft(Fz,Mz,3);                              % 1D fft in z
time = time + toc(ift_t);
x = x(1:Mx,1:My,1:Mz);                          % scale the result with numel

% switch back the dimensions
if per==1
    x = permute(x,[3 1 2]);
elseif per==2
    x = permute(x,[1 3 2]);
end

if(nargout==2)
    varargout{1} = time;
end

end
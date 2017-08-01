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


global_pad= setdiff(1:Mx,[local_pad zeromod]);                 % global pad vector

% periodic direction should be x!

time = 0;
ift_t = tic;
F0 = ifft2(x0mod,round(My*s0),round(Mz*s0)); % 2D fft on zero mod in y and z
Fr = ifft(ifft(xres,round(My*sl),2),round(Mz*sl),3); % 2D fft on
                                                     % local_pad in
                                                     % y and z
F = ifft(ifft(x,round(My*s),2),round(Mz*s),3);   % 2D fft on the rest on y and z
time = time + toc(ift_t);

if(~isempty(local_pad))
    Fxy(local_pad,:,:) = Fr(:,1:My,1:Mz);    % restrict to (local_pad,My,Mz)
end

Fxy(zeromod,:,:) = F0(1:My,1:Mz);           % restrict to
                                          % (zeromod,My,Mz)

Fxy(global_pad,:,:) = F(global_pad,1:My,1:Mz);  % update old values
                                                % and restrict

ift_t = tic;
x = ifft(Fxy,Mx);                              % 1D fft in x
time = time + toc(ift_t);
x = x(1:Mx,1:My,1:Mz);                          % scale the result with numel

if(nargout==2)
    varargout{1} = time;
end

end
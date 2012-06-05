function [phi A1 A2 A3] = stresslet_direct_real_fast( idx, x, f, nvec, xi, L, nbox, varargin)
% Ewald summation for the stresslet -- Real space part.
% Fast: saves interactions as matrix for subsequent iterations
%
% phi = stokes_ewald_direct_real( m, x, f, xi, op_A, nbox)
%        Evaluate potential, phi, at points x(idx). 
%        x    --  positions          N-by-3
%        nvec --  normal vectors   N-by-3
%        f    --  source strengths   N-by-3
%        xi   --  Ewald parameter
%        nbox --  periodic repications
%

VERBOSE = 0;

nosrc=size(f,1);
noeval=length(idx);  

tic
if nargin==10
    A1 = varargin{1};
    A2 = varargin{2};
    A3 = varargin{3};
else
    A1 = [];
    A2 = [];
    A3 = [];
end

if all( size(A1)==[noeval*3 nosrc] ) && ...
   all( size(A2)==[noeval*3 nosrc] ) && ...
   all( size(A3)==[noeval*3 nosrc] )
    cprintf(VERBOSE, '\tComputing real space sum using precomputed matrix.\n');
else

    p = -nbox:nbox; %boxes
    [p1 p2 p3] = ndgrid(p,p,p);
    p = [L(1)*p1(:) L(2)*p2(:) L(3)*p3(:)];
    Np = size(p,1);

    cprintf(VERBOSE, '\tComputing real space sum. Periodic images: %d\n', Np);

    % Keep interactions in three separate matrices, so that parfor can
    % slice them
    [A1 A2 A3] = deal(zeros(noeval*3,nosrc));
    idx = int32(idx);
    parfor n = 1:size(x,1) % particles
        
        xn = x(n,:); % source points
        nn = nvec(n,:);        
        
        % MEX inner loop
        tmp = stresslet_direct_real_mexcore(x,idx,xn,nn,n,nbox,xi,L);
        
        % MATLAB inner loop
%         tmp = zeros(noeval*3,3);
%         for m=1:noeval
%             xm = x(idx(m),:);
%             for j = 1:Np % periodic images
%                 if all(p(j,:)==0) && n==idx(m) % remove self interaction
%                   continue
%                 end
%                 a = stresslet_op_real( xm - xn + p(j,:), nn, xi);      
%                 k = m+noeval*[0 1 2];
%                 tmp(k,1:3) = tmp(k,1:3)+a;
%             end
%         end
        
        A1(:,n) = tmp(:,1);
        A2(:,n) = tmp(:,2);
        A3(:,n) = tmp(:,3);
    end
end

if VERBOSE
    fprintf('\t');
    toc
end

phi = reshape(A1*f(:,1) + A2*f(:,2) + A3*f(:,3), noeval, 3);

end

function [phi A] = stresslet_direct_real_fast( idx, x, f, nvec, xi, L, nbox, varargin)
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

VERBOSE = 1;

nosrc=size(f,1);
noeval=length(idx);  

tic
if nargin==8
    A = varargin{1};
else
    A = [];
end

if all( size(A)==[noeval*3 nosrc*3] )
    cprintf(VERBOSE, '\tComputing real space sum using precomputed matrix.\n');
else

    p = -nbox:nbox; %boxes
    [p1 p2 p3] = ndgrid(p,p,p);
    p = [L(1)*p1(:) L(2)*p2(:) L(3)*p3(:)];
    Np = size(p,1);

    cprintf(VERBOSE, '\tComputing real space sum. Periodic images: %d\n', Np);

    [A1 A2 A3] = deal(zeros(noeval*3,nosrc));
    parfor n = 1:size(x,1) % particles
        tmp = zeros(noeval*3,3);
        for j = 1:Np % periodic images
            for ii=1:noeval
                m=idx(ii); 
                if all(p(j,:)==0) && n==m % remove self interaction
                  continue
                end
                a = stresslet_op_real( x(m,:) - x(n,:) + p(j,:), nvec(n,:), xi);      
                k = ii+noeval*[0 1 2];
                tmp(k,1:3) = tmp(k,1:3)+a;
            end
        end
%       l = n+nosrc*[0 1 2];
%       A(:,l) = A(:,l)+tmp;
        A1(:,n) = tmp(:,1);
        A2(:,n) = tmp(:,2);
        A3(:,n) = tmp(:,3);
    end
    A = [A1 A2 A3];
end

if VERBOSE
    fprintf('\t');
    toc
end
phi = reshape(A*f(:), noeval, 3);

end

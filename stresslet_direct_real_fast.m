function [phi A1 A2 A3] = stresslet_direct_real_fast( idx, x, f, nvec, xi, L, nbox, rc, varargin)
% Ewald summation for the stresslet -- Real space part.
% Fast: saves interactions as matrix for subsequent iterations
%
% phi = stokes_ewald_direct_real( m, x, f, xi, op_A, nbox)
%        Evaluate potential, phi, at points x(idx). 
%        x    --  positions          N-by-3
%        nvec --  normal vectors   N-by-3
%        f    --  source strengths   N-by-3
%        xi   --  Ewald parameter
%        nbox --  periodic replications
%        rc   --  cutoff radius
%        [A1 A2 A3 sing_sub]

VERBOSE = 0;

nosrc=size(f,1);
noeval=length(idx);  

sing_sub = 0;

tic
if nargin>=11
    A1 = varargin{1};
    A2 = varargin{2};
    A3 = varargin{3};
    if nargin==12
        sing_sub = varargin{4};
    end
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
    cprintf(VERBOSE, '\tComputing real space sum matrix.\n');

    % Keep interactions in three separate matrices, one for each source
    % component
    idx = int32(idx);
    [A1 A2 A3] = stresslet_direct_real_mexcore(x,nvec,idx,nbox,rc,xi,L);
    
    if sing_sub
        cprintf(VERBOSE, 'Adding singularity subtraction to RS matrix\n');
        A1sum = sum(A1,2);
        A2sum = sum(A2,2);
        A3sum = sum(A3,2);
        for j=0:2
            for i=1:noeval
                rowno = i+noeval*j;
                A1(rowno,i) = A1(rowno,i)-A1sum(rowno);
                A2(rowno,i) = A2(rowno,i)-A2sum(rowno);
                A3(rowno,i) = A3(rowno,i)-A3sum(rowno);
            end
        end
    end
end

if VERBOSE
    fprintf('\t');
    toc
end

phi = reshape(A1*f(:,1) + A2*f(:,2) + A3*f(:,3), noeval, 3);

end

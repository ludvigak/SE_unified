function [phi A] = stresslet_direct_real_fast( idx, x, f, nvec, xi, L, nbox, rc, varargin)
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
%        [A sing_sub]

VERBOSE = 0;

nosrc=size(f,1);
noeval=length(idx);  

sing_sub = 0;

mytic=tic;
if nargin>=9
    A = varargin{1};
    if nargin==10
        sing_sub = varargin{2};
    end
else
    A = {[]};
end

% Just check size of one matrix
if all( size(A{1,1})==[noeval nosrc] )
    cprintf(VERBOSE, '\tComputing real space sum using precomputed matrix.\n');
else
    cprintf(VERBOSE, '\tComputing real space sum matrix.\n');

    % Keep interactions in three separate matrices, one for each source
    % component
    idx = int32(idx);
    A = stresslet_direct_real_mexcore(x,nvec,idx,nbox,rc,xi,L);
    
    if sing_sub
        cprintf(VERBOSE, 'Adding singularity subtraction to RS matrix\n');
        for k1=1:3
            for k2=k1:3
                Asum = sum(A{k1,k2},2);
                for i=1:noeval
                    A{k1,k2}(i,i) = A{k1,k2}(i,i) - Asum(i);
                end
            end
        end
        
    end
end

if VERBOSE
    fprintf('\t');
    toc(mytic)
end

phi = zeros(noeval,3);    
phi(:,1) = A{1,1}*f(:,1) + A{1,2}*f(:,2) + A{1,3}*f(:,3);
phi(:,2) = A{1,2}*f(:,1) + A{2,2}*f(:,2) + A{2,3}*f(:,3);
phi(:,3) = A{1,3}*f(:,1) + A{2,3}*f(:,2) + A{3,3}*f(:,3);

end

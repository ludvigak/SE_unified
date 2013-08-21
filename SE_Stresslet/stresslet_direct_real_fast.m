function [phi A] = stresslet_direct_real_fast( idx, x, f, nvec, xi, L, nbox, rc, varargin)
% Ewald summation for the stresslet -- Real space part.
% Fast: saves interactions as matrix for subsequent iterations
% Sums over layers
%
% [phi A] = stresslet_direct_real_fast( idx, x, f, nvec, xi, L, nbox, rc, 
%                                       [A, sing_sub, ns_quad_mtrx])
%        Evaluate potential, phi, at points x(idx). 
%        x    --  positions          N-by-3
%        f    --  source strengths   N-by-3
%        nvec --  normal vectors   N-by-3
%        xi   --  Ewald parameter
%        L    --  Box size [3]
%        nbox --  periodic replications 
%        rc   --  cutoff radius
%        A    --  Precomputed matrix
%        sing_sub -- Use singularity subtraction
%        ns_quad_mtrx -- Correction matrices for nearly singular quadrature
VERBOSE = 0;

nosrc=size(f,1);
noeval=length(idx);  
x = recenter_points(x, L);

sing_sub = 0;
ns_quad_mtrx = [];

A = {[]};

if VERBOSE
    mytic=tic;
end

nvarargin = numel(varargin);
if nvarargin
    if numel(varargin{1})
        A = varargin{1};
    end
    if nvarargin > 1
        sing_sub = varargin{2};
    end
    if nvarargin > 2
        ns_quad_mtrx = varargin{3};
    end
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
    
    % NOTE: No test suite for this in place yet
    if isstruct(ns_quad_mtrx)
        A = apply_ns_quad_mtrx(A, ns_quad_mtrx);
    end
    
    if sing_sub
        stic=tic;
        for k1=1:3
            for k2=k1:3
                Asum = sum(A{k1,k2},2);
                for i=1:noeval
                    A{k1,k2}(i,idx(i)) = A{k1,k2}(i,idx(i)) - Asum(i);
                end
            end
        end
            cprintf(VERBOSE, ...
                'Added singularity subtraction in %.3f seconds.\n',...
                toc(stic));
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

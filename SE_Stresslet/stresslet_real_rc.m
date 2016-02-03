function [res A varargout] = stresslet_real_rc( x, q, nvec, xi, box, rc, varargin)
% Stresslet real space summation with truncation radius rc
%
% Usage: 
%
% * [res]      = stresslet_real_rc( x, f, nvec, xi, box, rc);
% * [res AMAT] = stresslet_real_rc( x, f, nvec, xi, box, rc, [],   sing_sub, ns_quad_mtrx);
% * [res AMAT] = stresslet_real_rc( x, f, nvec, xi, box, rc, AMAT, sing_sub, ns_quad_mtrx);
% * [res AMAT R C V PER] = stresslet_real_rc( x, f, nvec, xi, box, rc);
%
% :param sing_sub: singularity subtraction, default 0. Only for matrix comp.
% :param rc: cutoff radius, must be <= min(box)/2
    
    VERBOSE = 0;
    
    check_inputs(box,rc);
    N = size(x,1);
    x = recenter_points(x, box);

    nvarargin = numel(varargin);
    if nvarargin>0
        A = varargin{1};
        if numel(A)==0 || numel(A{1})==0
            clear A;
        end
    end
    
    sing_sub = 0;
    if nvarargin>1
        sing_sub = varargin{2};
    end

    ns_quad_mtrx = [];
    if nvarargin > 2
        ns_quad_mtrx = varargin{3};
    end
    
    if ~exist('A','var')
        if nargout==1
            % Compute reults matrix-free
            if sing_sub
                error('Singularity subtraction not available for matrix-free call.')
            end
            if isstruct(ns_quad_mtrx)
                error('NS quad not available for matrix-free call.')
            end
            res = stresslet_real_rc_nomatrix_mex( x, nvec, q, box, rc, xi);
            return
        elseif nargout==2
            A = stresslet_real_rc_mex(x,nvec,box,rc,xi);
        elseif nargout==6
            [A R C V PER] = stresslet_real_rc_mex(x,nvec,box,rc,xi);
            [varargout{1:4}] = deal(R,C,V,PER);
        else
            error('Invalid usage');
        end
        
        % NOTE: No test suite for this in place yet
        if isstruct(ns_quad_mtrx)
            A = apply_ns_quad_mtrx(A, ns_quad_mtrx);
        end
        
        if sing_sub
            stic=tic;
            for k1=1:3
                for k2=k1:3
                    Asum = full(sum(A{k1,k2},2));
                    A{k1,k2} = A{k1,k2} - spdiags(Asum,0,N,N);
                end
            end
            cprintf(VERBOSE, ...
                '[RSRC] Added singularity subtraction in %.3f seconds.\n',...
                toc(stic));
        end
    end
    
    % Compute result using sparse matrices
    a = tic;
    res = zeros(N,3);    
    res(:,1) = A{1,1}*q(:,1) + A{1,2}*q(:,2) + A{1,3}*q(:,3);
    res(:,2) = A{1,2}*q(:,1) + A{2,2}*q(:,2) + A{2,3}*q(:,3);
    res(:,3) = A{1,3}*q(:,1) + A{2,3}*q(:,2) + A{3,3}*q(:,3);
    
    cprintf(VERBOSE,'[RSRC] matvec time %.3f seconds.\n',toc(a));
    
    if nargout==1
        A = [];
    end
    
    function check_inputs(box,rc)
        % Same algorithm as in C code, but with checks, to make sure boxing
        % is OK
        ncell = floor( min(box)/rc ); % number of cells in smallest direction
        rn = min(box) / ncell; % cell size
        ncell = box/rn; % number of cells
        % Make assertions
        errstr = sprintf(' (rc=%g, box=(%d,%d,%d), ncell=(%d,%d,%d))',...
                    rc, box(1),box(2),box(3),ncell(1),ncell(2),ncell(3));
        assert( all(ncell>1), ['Must be minimum 2 boxes in each direction' errstr]);
        assert( all(mod(ncell,1)==0 ), ['Box not divisible by cell size' errstr]);
    end
end
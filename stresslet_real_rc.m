function [res A varargout] = stresslet_real_rc( x, q, nvec, xi, box, rc, varargin)
    % Stresslet real space summation with truncation radius rc
    % Usage: 
    % [res]      = stresslet_real_rc( x, f, nvec, xi, box, rc);
    % [res AMAT] = stresslet_real_rc( x, f, nvec, xi, box, rc);
    % [res AMAT] = stresslet_real_rc( x, f, nvec, xi, box, rc, AMAT);
    % [res AMAT R C V PER] = stresslet_real_rc( x, f, nvec, xi, box, rc);
    
    VERBOSE = 1;
    
    check_inputs(box,rc);
    N = size(x,1);

    if nargin>6
        A = varargin{1};
    else
        if nargout==1
            % Compute reults matrix-free
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
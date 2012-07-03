function [res A varargout] = stresslet_real_rc( x, q, nvec, xi, box, rc, varargin)
    % Stresslet real space summation with truncation radius rc
    % Usage: 
    % [res AMAT] = stresslet_real_rc( x, f, nvec, xi, box, rc);
    % [res AMAT] = stresslet_real_rc( x, f, nvec, xi, box, rc, AMAT);
    % [res AMAT R C V PER] = stresslet_real_rc( x, f, nvec, xi, box, rc);
    
    check_inputs(box,rc);
    N = size(x,1);

    if nargin>6
        A = varargin{1};
    else
        if nargout==2
            A = stresslet_real_rc_mex(x,nvec,box,rc,xi);
        elseif nargout==6
            [A R C V PER] = stresslet_real_rc_mex(x,nvec,box,rc,xi);
            [varargout{1:4}] = deal(R,C,V,PER);
        else
            error('Invalid usage');
        end
    end

    res = zeros(N,3);    
    res(:,1) = A{1,1}*q(:,1) + A{1,2}*q(:,2) + A{1,3}*q(:,3);
    res(:,2) = A{1,2}*q(:,1) + A{2,2}*q(:,2) + A{2,3}*q(:,3);
    res(:,3) = A{1,3}*q(:,1) + A{2,3}*q(:,2) + A{3,3}*q(:,3);
    
    function check_inputs(box,rc)
        ncell = floor( min(box)/rc ); % number of cells in smallest direction
        rn = min(box) / ncell; % cell size
        ncell = box/rn; % number of cells
        % Make assertions
        assert( all(ncell>1), 'Must be minimum 2 boxes in each direction');
        assert( all(mod(ncell,1)==0 ), 'Box not divisible by cell size' );
    end
end
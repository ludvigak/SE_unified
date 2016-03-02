function [ xvec ] = recenter_points( xvec, box )
% [ xvec ] = recenter_points( xvec, box )
%
% Recenter points in mother box

    if numel(xvec)==0
        return;
    end
    
    xvec = bsxfun(@plus, xvec, box);
    xvec = bsxfun(@mod, xvec, box);

end


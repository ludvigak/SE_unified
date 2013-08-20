function [ xvec ] = recenter_points( xvec, box )
% [ xvec ] = recenter_points( xvec, box )
%
% Recenter points in mother box

xvec = bsxfun(@mod, xvec, box);

end


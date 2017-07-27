function [x, f] = vector_system(N, box, varargin)
%VECTOR_SYSTEM vector systems of different dimensions.
%   VECTOR_SYSTEM(N, [L1 L2 L3])=VECTOR_SYSTEM(N, [L1 L2 L3],1) generates
%   a 1d random system on a box of size [L1 L2 L3]. 
%
%   VECTOR_SYSTEM(N, [L1 L2 L3],DIM) generates a DIM dimensional random system.

if isempty(varargin)
    dim = 1;
else
    dim = varargin{1};
end

x = bsxfun(@times, rand(N, 3), box);
f = 1-2*rand(N, dim);

if dim==1
    f = f - repmat( mean(f,1), N, 1); % neutrality
end
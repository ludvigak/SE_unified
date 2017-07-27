function [x, f] = vector_system(N, box, varargin)

if isempty(varargin)
    dim = 3;
else
    dim = varargin{1};
end

x = bsxfun(@times, rand(N, 3), box);
f = 1-2*rand(N, dim);

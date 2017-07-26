function [x, f] = laplace_system(N, box)
x = bsxfun(@times, rand(N, 3), box);
%f = 1-2*rand(N, 1);
f = rand(N,1);
f = f - mean(f);


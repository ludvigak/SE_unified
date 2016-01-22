function [x, t, xe] = generate_state(N, box)
x = bsxfun(@times, box, rand(N,3));
t = 1-2*rand(N,3);
xe = bsxfun(@times, box, rand(N,3));

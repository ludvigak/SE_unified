function u = rotlet(x, y, t)
% u = rotlet(x, y, t(y))
%
% Compute rotlet potential
% u_j = (t x r)_i / r^3

r = bsxfun(@minus, x, y);
ri3 = sum(r.^2, 2).^(-3/2);

if size(y, 1) == 1
    % Vectorize in x
    txr = bsxfun(@cross, t, r, 1);
    u = bsxfun(@times, txr, ri3);
else
    error('y input must be point');
end

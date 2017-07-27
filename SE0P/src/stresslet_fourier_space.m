function [u varargout] = stresslet_fourier_space(x, f, opt, pre)

q = f(:,[1 2 3]);
n = f(:,[4 5 6]);

N = size(q, 1);
s = zeros(N, 3, 3);
for j=1:3
    for l=1:3
        s(:, j, l) = q(:, j).*n(:, l);
    end
end

if nargout == 2
   [u walltime] = fse_fourier_space(x, s, opt, pre, @stresslet_k_scaling);
   varargout{1} = walltime;
else
   u = fse_fourier_space(x, s, opt, pre, @stresslet_k_scaling);
end
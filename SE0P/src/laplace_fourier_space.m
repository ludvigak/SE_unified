function u = laplace_fourier_space(x, f, opt, pre,varargin)

if(nargin==4)
    u = fse_fourier_space(x, f, opt, pre, @laplace_k_scaling);
else
    u = fse_fourier_space(x, f, opt, pre, @ ...
                          laplace_k_scaling_force);
end

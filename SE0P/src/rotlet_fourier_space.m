function [u varargout]= rotlet_fourier_space(x, f, opt, pre)

if nargout==2
   [u walltime]= fse_fourier_space(x, f, opt, pre, @rotlet_k_scaling);
    varargout{1} = walltime;
else
   u = fse_fourier_space(x, f, opt, pre, @rotlet_k_scaling);
end
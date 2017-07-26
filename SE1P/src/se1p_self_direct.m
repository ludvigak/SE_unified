function [u varargout]= se1p_self_direct(idx,x,f,opt,varargin)
    N = numel(f);
    xi = opt.xi;
if nargin ==5
    MATLAB = varargin{1};
else
    MATLAB = false;
end    

time = tic;
if(MATLAB)
    c = 2*opt.xi/sqrt(pi);
    u = zeros(numel(idx),1);
    u(idx) = -c*f(idx);
else
    u = SE1P_direct_self_mex(idx,x,f,opt);
end
walltime = toc(time);

if(nargout==2)
    varargout{1} = walltime;
end 

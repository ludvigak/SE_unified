function [u varargout] = rotlet_real_space(x, f, opt)

xi = opt.xi;
rc = opt.rc;

N = size(x, 1);
timestamp = tic();
[idx, d] = rangesearch(x, x, rc);
walltime.nblist = toc(timestamp);

timestamp = tic();
MATLAB = false;
if MATLAB
    u = zeros(N, 3);
    for target = 1:N
        nnb_idx = idx{target}(2:end);
        rvec = bsxfun(@minus, x(target,:), x(nnb_idx,:));
        r = d{target}(2:end)';
        fsrc = f(nnb_idx,:);
        
        r2 = r.^2;
        fxr = cross(fsrc, rvec, 2);
        
        A = ( erfc(xi*r)./r + 2*xi*exp(-xi^2*r2)/sqrt(pi) )./r2;
        u(target,:) = sum( bsxfun(@times, fxr, A), 1);
    end
else
    u = SE0P_rotlet_rsrc_mex(x,f,idx,d,xi);
end
walltime.eval = toc(timestamp);

if nargout==2
    walltime.total = walltime.nblist + walltime.eval;
    varargout{1} = walltime;
end

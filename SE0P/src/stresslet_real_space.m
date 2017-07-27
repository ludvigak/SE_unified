function [u varargout] = stresslet_real_space(x, f, opt)

q = f(:,[1 2 3]);
n = f(:,[4 5 6]);

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
        qsrc = q(nnb_idx,:);
        nsrc = n(nnb_idx,:);
        
        r2 = r.^2;
        c = xi^2*r2;    
        % Beenakker
        %C=-2./r2.^2.*( 3.0./r.*erfc(xi*r) + 2.0*xi/sqrt(pi)*(3.0+2.0*c-4.0*c.^2).*exp(-c) );
        %D=8/sqrt(pi)*xi^3*(2.0-c).*exp(-c);
        % Hasimoto
        C=-2./r2.^2.*( 3.0./r.*erfc(xi*r) + 2.0*xi/sqrt(pi)*(3.0+2.0*c).*exp(-c) );
        D=4/sqrt(pi)*xi^3.*exp(-c);
        
        rdotn = sum(rvec .* nsrc, 2);
        rdotq = sum(rvec .* qsrc, 2);
        ndotq = sum(nsrc .* qsrc, 2);
        
        u(target,:) = sum(...
            bsxfun(@times, C.*rdotn.*rdotq + D.*ndotq, rvec) + ...
            bsxfun(@times, D.*rdotq, nsrc) +...
            bsxfun(@times, D.*rdotn, qsrc) ...
            , 1);
    end
else
    u = SE0P_stresslet_rsrc_mex(x,f,idx,d,xi);
end
walltime.eval = toc(timestamp);

if nargout==2
    walltime.total = walltime.nblist + walltime.eval;
    varargout{1} = walltime;
end

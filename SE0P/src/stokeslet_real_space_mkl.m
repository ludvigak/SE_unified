function [u varargout] = stokeslet_real_space_mkl(x, f, opt)

walltime.nblist = 0;

timestamp = tic();
u = SE0P_stokeslet_rsrc_cell_mkl_mex(x', f', opt.rc, opt.xi, opt.box)';
walltime.eval = toc(timestamp);

if nargout==2
    walltime.total = walltime.nblist + walltime.eval;
    varargout{1} = walltime;
end

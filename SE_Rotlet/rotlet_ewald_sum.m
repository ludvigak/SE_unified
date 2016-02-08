function u = rotlet_ewald_sum(xe, x, t, opt)
% Ewald summation for 3P rotlet.
%
% :param xe: target locations
% :param x: source locaitons
% :param t: source strengths
% :param opt: Ewald options
% :param opt.box: box size
% :param opt.xi: Ewald parameter :math:`\xi`
% :param opt.M: Fourier space grid size
% :param opt.P: Gaussian width
% :param opt.rc: Real space cutoff radius


uf = SE_Rotlet(xe, x, t, opt.xi, opt);
ur = rotlet_direct_rsrc(xe, x, t, opt);

u = uf + ur;
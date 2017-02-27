function u = rotlet_direct_rsrc(xe, x, t, opt)
% Rotlet real space direct summation with cutoff radius, MEX implementation.
%
% :param xe: target locations
% :param x: source locaitons
% :param t: source strengths
% :param opt: Ewald options
% :param opt.box: box size
% :param opt.xi: Ewald parameter :math:`\xi`
% :param opt.rc: Cutoff radius


assert(opt.rc <= min(opt.box), 'Required: rc <= min(box)')
assert(opt.xi > 0, 'Required: xi > 0')
assert(size(xe, 2) == 3, 'xe must be Mx3');
assert(size(x, 2) == 3, 'x must be Nx3');
assert(all(size(x) == size(t)), 'x and t must have same dimensions');

u = rotlet_direct_rsrc_mex(xe, x, t, opt);


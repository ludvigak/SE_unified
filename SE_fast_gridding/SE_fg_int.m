function u = SE_fg_int(x, U, opt)
% u = SE_fg_int(x, U, opt)
%   Return integration at x of on-grid values U
%
% u = SE_fg_int(x, {U1, ..., UN}, opt)
%   Return integration at x of on-grid values U1,...,UN
%   u has N columns
%

% parameters and constants
opt = parse_params(opt);

% Integrator
SI = SE_FGG_precomp(x, opt.xi, opt);
iperm = @(u) u(SI.iperm,:);
int_fcn = @(F) iperm(SE_fg_int_split_mex(0,F,opt,SI.zs,SI.zx,SI.zy,SI.zz,SI.idx));

if iscell(U)
    u = zeros(size(x, 1), numel(U));
    for i=1:numel(U)
        u(:, i) = int_fcn(U{i});
    end
else
    u = int_fcn(U);
end


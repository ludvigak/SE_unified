function H = SE_fg_extend_fcn(F, opt)

% parameters and constants
opt = parse_params(opt);

% Call MEX
if iscell(F)
    H = cell(size(F));
    for i=1:numel(F)
        H{i} = SE_fg_extend_fcn_mex(F{i}, opt); 
    end
else
    H = SE_fg_extend_fcn_mex(F, opt);
end

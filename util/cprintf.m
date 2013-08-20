function cprintf(verb, fmt, varargin)
% Conditional fprintf
% cprintf(true, fmt, args) is equivalent to fprintf(fmt, args)
% otherwise silent. Typical usage:
% 
% cprintf(verb, 'status message')
%  
% instead of 
%
% if verb, fprintf('status message'), end
%
% Dag Lindbo, dag@kth.se, 2009-02-11

% Print nothing
if ~verb, return, end

% pass args to printf
fprintf(fmt, varargin{:});

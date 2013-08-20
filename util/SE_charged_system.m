function [x h] = SE_charged_system(N,L,s)
% Charge neutral system och particles, 'scalar' 
%   or 'vector' valued charges.

x = repmat(L,N,1).*rand(N,3);
switch s
  case 'scalar' 
    % Laplace
    h = rand( N, 1);
    h = h - mean(h); % neutrality
  case 'vector' 
    % Stokes
    h = rand(N,3);
    h = h - repmat( mean(h,1), N, 1); % neutrality
end
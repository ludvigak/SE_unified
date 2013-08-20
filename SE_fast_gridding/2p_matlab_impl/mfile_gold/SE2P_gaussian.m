function [Q ix iy iz] = SE2P_gaussian(x,opt,c)

h = opt.h;
P = opt.P;
M = opt.M;
%a = opt.wbox(3,1);
a = opt.a;

if isodd(P)
    [ix px] = odd_grid_restrict_0_wrap(x(1),h,P,M);
    [iy py] = odd_grid_restrict_0_wrap(x(2),h,P,M);
    % z: on staggered grid, starting at a.
    [iz pz] = odd_grid_restrict_a_stagger(x(3),h,P,a);
else
    [ix px] = even_grid_restrict_0_wrap(x(1),h,P,M);
    [iy py] = even_grid_restrict_0_wrap(x(2),h,P,M);
    % z: on staggered grid, starting at a.
    [iz pz] = even_grid_restrict_a_stagger(x(3),h,P,a);
end

% distances
[DX DY DZ] = ndgrid(px-x(1),py-x(2),pz-x(3));
D2 = DX.^2 + DY.^2 + DZ.^2;

% gaussian
Q = (c/pi)^1.5*exp(-c*D2);

function [ix px] = even_grid_restrict_0_wrap(x,h,P,M)

idx = floor(x/h)+1;
P_half = P/2;
ix = (idx-(P_half-1)):(idx+P_half);
px = h*(ix-1);
ix = mod(ix-1,M)+1;

function [ix px] = odd_grid_restrict_0_wrap(x,h,P,M)

idx = round(x/h)+1;
P_half = (P-1)/2;
ix = (idx-P_half):(idx+P_half);
px = h*(ix-1);
ix = mod(ix-1,M)+1;

function [ix px] = even_grid_restrict_a_stagger(x,h,P,a)

idx = floor((x-(a+h/2))/h)+1;
P_half = P/2;
ix = (idx-(P_half-1)):(idx+P_half);
px = h*(ix-1) + a + h/2;

function [ix px] = odd_grid_restrict_a_stagger(x,h,P,a)

idx = round((x-(a+h/2))/h)+1;
P_half = (P-1)/2;
ix = (idx-P_half):(idx+P_half);
px = h*(ix-1) + a + h/2;

function q = isodd(P)
q=mod(P,2)==1;
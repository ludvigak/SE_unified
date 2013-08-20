function [Q ix iy iz] = SE_gaussian(x,h,M,P,c)

% support indices
idx = round(x/h)+1;
idx_low = idx-(P-1)/2;
idx_upp = idx+(P-1)/2;

supp_x = idx_low(1):idx_upp(1);
supp_y = idx_low(2):idx_upp(2);
supp_z = idx_low(3):idx_upp(3);

[SX SY SZ] = ndgrid(supp_x,supp_y,supp_z);

% grid points (assume interval start at 0)
PX = h*(SX-1);
PY = h*(SY-1);
PZ = h*(SZ-1);

% distances
D2 = (PX-x(1)).^2 + (PY-x(2)).^2 + (PZ-x(3)).^2;

% assignment indices (wrapped to produce periodicity)
ix = mod(supp_x-1,M(1))+1;
iy = mod(supp_y-1,M(2))+1;
iz = mod(supp_z-1,M(3))+1;

% gaussian
Q = (c/pi)^1.5*exp(-c*D2);

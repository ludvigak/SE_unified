clear

% Test for small N to include MATLAB code in comparison
fgg_timing(11, 1000)

 % Test for some varying P to make sure SIMD codes are working
fgg_timing([12 13 14 15 16], 2000)

% Finish off with a large N to show difference between methods
fgg_timing(16, 1e5)
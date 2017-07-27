function GR = biharmonic(K, R)
% Fourier transform of r truncated at R

K2 = K.^2;
GR = 8*pi*( ...
    (2-R^2*K2).*cos(R*K) + 2*R*K.*sin(R*K) - 2 ...
    ) ./ (2*K2.^2);
GR(K==0) = pi*R^4;

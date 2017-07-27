function GR = harmonic(K, R)
% Fourier transform of 1/r truncated at R

GR = 8*pi*(sin(K*R/2)./K).^2;
GR(K==0) = 2*pi*R^2;

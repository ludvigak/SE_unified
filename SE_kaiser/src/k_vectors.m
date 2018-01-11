function [k1 k2 k3] = k_vectors(M,box)
if (all(mod(M,2)==0))
  MM = M/2;
  k1 = (2*pi/box(1))*[0:(MM(1)-1) -MM(1):-1];
  k2 = (2*pi/box(2))*[0:(MM(2)-1) -MM(2):-1];
  k3 = (2*pi/box(3))*[0:(MM(3)-1) -MM(3):-1];

elseif(all(mod(M-1,2)==0))
  MM = (M-1)/2;
  k1 = (2*pi/box(1))*[0:MM(1) -MM(1):-1];
  k2 = (2*pi/box(2))*[0:MM(2) -MM(2):-1];
  k3 = (2*pi/box(3))*[0:MM(3) -MM(3):-1];

else error('k-vectors not computed (FIXME)');
end
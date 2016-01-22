function E = rotlet_est_real_long(t, box, rc, xi)

Q = sum(t(:).^2);
V = prod(box);

base = erfc(xi*rc).^2 + sqrt(2/pi)*rc*xi.*erfc(sqrt(2)*rc*xi);
E = sqrt(8*pi*Q./(3*V*rc).*base);

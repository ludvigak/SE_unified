function E = rotlet_est_real(t, box, rc, xi)

Q = sum(t(:).^2);
V = prod(box);
E = sqrt(8*Q./(3*V*rc)).*exp(-(xi*rc).^2);
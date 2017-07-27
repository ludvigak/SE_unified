function rc=calc_rc(xi, tol)

rc = 5;
i = 0;

while (erfc(rc*xi) > tol);
    i = i + 1;
    rc = rc*2;
end

n    = i+60; % search tolerance is 2^-60
low  = 0;
high = rc;
for i = 0:n-1
    rc = (low+high)/2;
    if (erfc(rc*xi) > tol)
        low = rc;
    else
        high = rc;
    end
end

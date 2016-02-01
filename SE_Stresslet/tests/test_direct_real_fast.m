clear

format long e
rand('state',1)

fprintf('\n\n')

box = [1.1 1+pi/10 sqrt(2)/2];
xi = 1;

NOL = 2; % layers
rc = sqrt(box*box')*(NOL+1);
N = 15;

[x f nvec] = generate_state(N,box);
idx = 1:N;

fprintf('Real space summation of %d points over %d layers.\n',N,NOL);

tic;
fprintf('MATLAB direct summation: \n\t')
uref = stresslet_direct_real(idx, x, f, nvec, xi,  box, NOL);
toc
wref = toc;
tic
fprintf('MEX matrix assembly + matvec: \n\t')
[ufst A] = stresslet_direct_real_fast(idx, x, f, nvec, xi,  box, NOL, rc);
toc
wfst = toc;
timings = [wref wfst]

res1 = abs(ufst-stresslet_direct_real_fast(idx, x, f, nvec, xi,  box, NOL, rc, A));
if max(res1(:))>1e-10
    error('EWALD FAST RS: FAILED MATRIX')
end


ufst2 = zeros(N,3);
for i=idx
    ufst2(i,:) = stresslet_direct_real_fast(i, x, f, nvec, xi,  box, NOL, rc);
end
res3 = ufst(:)-ufst2(:);
if norm(res3,inf) > 1e-13
    error('EWALD FAST RS: PARTIAL SOURCES')
end

res2 = abs(uref-ufst);
fprintf('Max diff: %g\n', max(res2(:)));

if max(res2(:))>1e-10
    maxres2 = max(res2(:))
    error('EWALD FAST RS: FAILED')
else
    fprintf('\n********** EWALD FAST RS: OK **********\n\n')
end


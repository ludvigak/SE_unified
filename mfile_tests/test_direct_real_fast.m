format long e
rand('state',1)

if matlabpool('size')==0
    matlabpool('open')
end

fprintf('\n\n')

box = [1.1 1+pi/10 sqrt(2)/2];
xi = 1;

NOL = 2; % layers
N = 15;

[x f nvec] = generate_state(N,box);
idx = 1:N;

fprintf('Real space summation of %d points over %d layers.\n',N,NOL);

tic;
% fprintf('MATLAB direct summation: \n\t')
% uref = stresslet_direct_real(idx, x, f, nvec, xi,  box, NOL);
% toc
wref = toc;
tic
fprintf('MEX matrix assembly + matvec: \n\t')
[ufst A] = stresslet_direct_real_fast(idx, x, f, nvec, xi,  box, NOL);
toc
wfst = toc;
timings = [wref wfst]

res1 = abs(ufst-stresslet_direct_real_fast(idx, x, f, nvec, xi,  box, NOL, A));
if max(res1(:))>1e-10
    error('EWALD FAST RS: FAILED MATRIX')
end


res2 = abs(uref-ufst);
if max(res2(:))>1e-10
    error('EWALD FAST RS: FAILED')
end
fprintf('Max diff: %g\n', max(res2(:)));

fprintf('\n********** EWALD FAST RS: OK **********\n\n')


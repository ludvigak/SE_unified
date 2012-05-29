format long e
rand('state',1)

if matlabpool('size')==0
    matlabpool('open')
end

box = [1.1 1+pi/10 sqrt(2)/2];
xi = 1;

NOL = 3; % layers
N = [400];
timings = [];
for i=1:length(N)
    [x f nvec] = generate_state(N(i),box);
    idx = 1:N(i);
    tic;
    uref = 0;
%     uref = stresslet_direct_real(idx, x, f, nvec, xi,  box, NOL);
    wref = toc;
    tic
    [ufst A] = stresslet_direct_real_fast(idx, x, f, nvec, xi,  box, NOL);
    wfst = toc;
    timings(end+1,1:2) = [wref wfst];
end

timings

res1 = abs(ufst-stresslet_direct_real_fast(idx, x, f, nvec, xi,  box, NOL, A));
if max(res1(:))>1e-10
    error('EWALD FAST RS: FAILED MATRIX')
end


res2 = abs(uref-ufst);
if max(res2(:))>1e-10
    error('EWALD FAST RS: FAILED')
end

fprintf('\n********** EWALD FAST RS: OK **********\n\n')


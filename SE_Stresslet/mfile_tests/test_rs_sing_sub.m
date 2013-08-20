function status = test_rs_sing_sub()

N = 1000;
box = [1 1 1];
rc = 0.5;
xi = 5;
sing_sub = 1;

[x f n] = generate_state(N,box);
idx = 5:3:N-13;

% rs rc
[res_rsrc ~] = stresslet_real_rc( x, f, n, xi, box, rc, [], sing_sub);
res_rsrc = res_rsrc(idx,:);
% direct w. matrix
res_direct = stresslet_direct_real_fast(idx, x, f, n, xi,  box, 1, rc, [], sing_sub);
% brute
res_brute = zeros(numel(idx),3);
parfor j = 1:numel(idx)
    i = idx(j);
    fvec_sub = bsxfun(@minus,f,f(i,:));
    res_brute(j,:) = ...
        stresslet_direct_real_fast(i,x,fvec_sub,n,...
        xi,box,1,rc,{[]},0);
end

e1 = res_brute - res_rsrc;
e1_inf = norm(e1(:),inf) / norm(res_brute(:),inf);

e2 = res_brute - res_direct;
e2_inf = norm(e2(:),inf) / norm(res_brute(:),inf);

fprintf('e_rsrc = %.2e\ne_direct = %.2e\n',e1_inf,e2_inf);

e_inf = max(e1_inf,e2_inf);

if e_inf>1e-12
    warning('Test:RSSingularitySubtraction','Real space singularity subtraction failed.');
    status = 0;
else
    disp('Real space singularity subtraction OK');
    status = 1;
end
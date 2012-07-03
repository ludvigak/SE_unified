% Test completely MEXed real space rc based on cell lists

disp('==== test_real_rc ====')

addpath('~/workspace/stresslet_ewald/matlab')

N = 3000;
box = [3 2 1.5];
xi = 1;
rc = 0.4;
NOL = 1;

[x f nvec] = generate_state(N,box);
idx = 1:N;

fprintf('* Direct rsrc\n\t')
a=tic;
[ufst A1 A2 A3] = stresslet_direct_real_fast(idx, x, f, nvec, xi,  box, NOL, rc);
toc(a)

ic = 1;
j = find(A1(1:N,ic));

figure(10)
drawBox(box,[0 0 0]), hold on
plot3(x(ic,1),x(ic,2),x(ic,3),'or');
plot3(x(j,1),x(j,2),x(j,3),'og');
drawnow


disp('* Real rc')
a=tic;
[res AMAT R C V PER] = stresslet_real_rc( x, f, nvec, xi, box, rc);
fprintf('\t')
toc(a)

j = find(AMAT{1,1}(ic,:));
plot3(x(j,1),x(j,2),x(j,3),'.b');
plot3(x(ic,1),x(ic,2),x(ic,3),'.b');
drawnow


fprintf('Matlab quicksort\n\t')
R = double(R);
C = double(C);
a=tic;
[~, perm] = sortrows([R C],[2 1]);
toc(a)

R = R(PER);
C = C(PER);
Ac = cell(3,3);
fprintf('Matlab assembly\n\t')
a=tic;
for i=1:3
    for j=i:3
        Ac{i,j} = sparse(R,C,V{i,j}(PER));
    end
end
toc(a)

% TESTS

% Check permutations
assert(all(perm==PER),'Permutations failed!')
disp('Permutations OK')

% Check matrices equal
Aref = sparse([A1 A2 A3]);
Anew = ...
    [Ac{1,1} Ac{1,2} Ac{1,3}
     Ac{1,2} Ac{2,2} Ac{2,3}
     Ac{1,3} Ac{2,3} Ac{3,3}];
relerr = norm(Aref-Anew,inf)/norm(Aref,inf) ;
assert(relerr<1e-14,'Explicit matrix differes from reference matrix!')

Anew = ...
    [AMAT{1,1} AMAT{1,2} AMAT{1,3}
     AMAT{1,2} AMAT{2,2} AMAT{2,3}
     AMAT{1,3} AMAT{2,3} AMAT{3,3}];
relerr2 = norm(Aref-Anew,inf)/norm(Aref,inf) ;
assert(relerr2<1e-14,'Returned matrix differs from reference matrix!')

disp('Matrices OK')

% Check results equal
uref = stresslet_direct_real_fast(idx, x, f, nvec, xi,  box, NOL, rc,A1,A2,A3);
u = stresslet_real_rc( x, f, nvec, xi, box, rc, AMAT);
udiff = norm(u-uref,inf)/norm(uref,inf);
assert(udiff<1e-14,'Results from RS methods differ!')

disp('Output OK');

clear Aref Anew Ac AMAT A1 A2 A3 R C V PER

disp('======================')
disp('********** EWALD RS RC: OK **********')

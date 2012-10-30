% Test completely MEXed real space rc based on cell lists


disp('==== test_real_rc ====')

addpath('~/workspace/stresslet_ewald/matlab')

N = 1500;
box = [3 2 1];
xi = 1;
rc = 0.4;
NOL = 1;
ic = 1;

rc = 0.4;


[x f nvec] = generate_state(N,box);
idx = 1:N;

disp('* Real rc')
a=tic;
[res AMAT R C V PER] = stresslet_real_rc( x, f, nvec, xi, box, rc);
fprintf('\t')
toc(a)

sfigure(9); clf
plot3(x(:,1),x(:,2),x(:,3),'.b'), hold on
plot3(x(ic,1),x(ic,2),x(ic,3),'or');
drawBox(box,[0 0 0]), hold on
drawnow

sfigure(10); clf
j = find(AMAT{1,1}(ic,:));
plot3(x(j,1),x(j,2),x(j,3),'.b'), hold on
plot3(x(ic,1),x(ic,2),x(ic,3),'.b');
drawBox(box,[0 0 0]), hold on
view(3)
drawnow

j1 = j;

Ac = cell(3,3);
fprintf('Matlab assembly\n\t')
a=tic;
R = double(R);
C = double(C);
for i=1:3
    for j=i:3
        Ac{i,j} = sparse(R,C,V{i,j}(PER),N,N);
    end
end
toc(a)

fprintf('* Direct rsrc\n\t')
a=tic;
[ufst A] = stresslet_direct_real_fast(idx, x, f, nvec, xi,  box, NOL, rc);
toc(a)

A1N = sparse(A{1,1});

j = find(A{1,1}(:,ic));
j2 = j;

drawBox(box,[0 0 0]), hold on
plot3(x(ic,1),x(ic,2),x(ic,3),'or');
plot3(x(j,1),x(j,2),x(j,3),'og');
drawnow

% TESTS

% Check sparsity pattern
assert( all( find(A{1,1})==find(AMAT{1,1}) ), 'Sparsity pattern differs!');

% Check matrices equal
Aref = sparse(...
    [A{1,1} A{1,2} A{1,3}
     A{1,2} A{2,2} A{2,3}
     A{1,3} A{2,3} A{3,3}]);
Anew = ...
    [Ac{1,1} Ac{1,2} Ac{1,3}
     Ac{1,2} Ac{2,2} Ac{2,3}
     Ac{1,3} Ac{2,3} Ac{3,3}];
relerr = norm(Aref-Anew,inf)/norm(Aref,inf) ;
assert(relerr<1e-14,'Explicit matrix differs from reference matrix!')

Anew = ...
    [AMAT{1,1} AMAT{1,2} AMAT{1,3}
     AMAT{1,2} AMAT{2,2} AMAT{2,3}
     AMAT{1,3} AMAT{2,3} AMAT{3,3}];
relerr2 = norm(Aref-Anew,inf)/norm(Aref,inf) ;
assert(relerr2<1e-14,'Returned matrix differs from reference matrix!')

disp('Matrices OK')

% Check results equal
uref = stresslet_direct_real_fast(idx, x, f, nvec, xi,  box, NOL, rc, A);
u = stresslet_real_rc( x, f, nvec, xi, box, rc, AMAT);
udiff = norm(u-uref,inf)/norm(uref,inf);
assert(udiff<1e-14,'Results from RS methods differ!')

disp('Output OK');

clear Aref Anew Ac AMAT A1 A2 A3 R C V PER

disp('======================')
disp('********** EWALD RS RC: OK **********')

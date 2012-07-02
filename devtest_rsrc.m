clear

% TODO: Save this as a complete test

rand('state',2)
addpath('~/workspace/stresslet_ewald/matlab')

N = 3000;
box = [1 1 1];
xi = 1;
rc = 0.3;
NOL = 1;

[x f nvec] = generate_state(N,box);
idx = 1:N;

tic
[ufst A1 A2 A3] = stresslet_direct_real_fast(idx, x, f, nvec, xi,  box, NOL, rc);
toc

A11 = A1(1:N,1:N);

clf
drawBox(box,[0 0 0]), hold on
% plot3(x(:,1),x(:,2),x(:,3),'.'), hold on

ic = 1;
plot3(x(ic,1),x(ic,2),x(ic,3),'or');

j = find(A11(:,ic));
plot3(x(j,1),x(j,2),x(j,3),'og');

% % Check (only small rc)
xr = [x(:,1)-x(ic,1) x(:,2)-x(ic,2) x(:,3)-x(ic,3)];
r = sqrt(sum(xr.^2,2));
j2 = find(r<rc);
j2 = j2(j2~=ic) ;

%%

a=tic;
[AMAT R C V PER] = stresslet_real_rc(x,nvec,box,rc,xi);
mytime = toc(a)

R = double(R);
C = double(C);

a=tic;
[~, perm] = sortrows([R C],[2 1]);
qsorttime = toc(a)

assert(all(perm==PER),'Permutations failed')



disp('Assembly sparse2')
tic
Ac = cell(3,3);
for i=1:3
    for j=i:3
        Ac{i,j} = sparse2(R,C,V{i,j});
    end
end
toc

disp('Assembly')
tic
R = R(PER);
C = C(PER);
Ac = cell(3,3);
for i=1:3
    for j=i:3
        Ac{i,j} = sparse(R,C,V{i,j}(PER));
    end
end
toc

j = find(Ac{1,1}(ic,:));
plot3(x(j,1),x(j,2),x(j,3),'.b');
plot3(x(ic,1),x(ic,2),x(ic,3),'.b');


% Check
Aref = sparse([A1 A2 A3]);
Anew = ...
    [Ac{1,1} Ac{1,2} Ac{1,3}
     Ac{1,2} Ac{2,2} Ac{2,3}
     Ac{1,3} Ac{2,3} Ac{3,3}];
 
relerr = norm(Aref-Anew,inf)/norm(Aref,inf) 
assert(relerr<1e-14)

% Check 2
Anew = ...
    [AMAT{1,1} AMAT{1,2} AMAT{1,3}
     AMAT{1,2} AMAT{2,2} AMAT{2,3}
     AMAT{1,3} AMAT{2,3} AMAT{3,3}];
 
relerr2 = norm(Aref-Anew,inf)/norm(Aref,inf) 
assert(relerr2<1e-14)


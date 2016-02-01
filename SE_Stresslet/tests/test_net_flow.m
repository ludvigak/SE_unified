clear
% Test that computes flow through bottom xy-plane is equal to analytical
% result.

format long
box = [3 1 2];

% Random setup
rng('default')
rng('shuffle')
xvec = [rand(1,2) 0.3+0.4*rand()].*box
fvec = rand(1,3)
nvec = rand(1,3)

% Setup params
TOL_r=1e-200;
NOL_r = 1;
xi = 20;
SE_opt.xi = xi;
SE_opt.box = box;
SE_opt.M = box*70;
SE_opt.P = 16;

% Discretize xy plan
N = 20;
Ngrid = box(1:2)*N;
h = 1/N;
xp = linspace(0, box(1)-h, Ngrid(1));
yp = linspace(0, box(2)-h, Ngrid(2));
zp = rand()*0.1;
hx = xp(2)-xp(1);
hy = yp(2)-yp(1);
[X Y] = meshgrid(xp,yp);
Z = zeros(size(X))+zp;

% Setup quad weights, quadrature rule works great because we are
% integrating periodic part.
Q = ones(size(X))*hx*hy;
assert(abs(sum(Q(:))-box(1)*box(2))<1e-10);

% Ewald summation
eval_x = [X(:) Y(:) Z(:)];
SE_opt.eval_x = eval_x;

SE_static  = SE_Stresslet_pre(xvec,xi,SE_opt);
[u_f stats] = SE_Stresslet([],xvec,fvec,nvec,xi,SE_opt, SE_static);

u_r = stresslet_direct_real([], xvec, fvec, nvec, xi,  box, NOL_r, TOL_r, eval_x);

xvec = [xvec; eval_x];
fvec = [fvec; eval_x*0];
nvec = [nvec; eval_x*0];
idx = 2:numel(X)+1;

u_z=stresslet_direct_fd_zero(idx, xvec, fvec, nvec, box);

%
u = u_r + u_f + u_z;

W_r = reshape(u_r(:,3),size(X));
I_r = sum(W_r(:).*Q(:));

W_f = reshape(u_f(:,3),size(X));
I_f = sum(W_f(:).*Q(:));

W_z = reshape(u_z(:,3),size(X));
I_z = sum(W_z(:).*Q(:));

W = reshape(u(:,3),size(X));
I = sum(W(:).*Q(:));

% Attempt I_f computation
zpoint = xvec(1,3);
zplane = zp;
rz = zpoint-zplane;

% Analytical formuls
U_f = 8*pi/box(3)*dot(nvec(1,:),fvec(1,:))*(box(3)/2-mod(rz,box(3)));

disp(' ')
disp('I_f/U_f')
disp([I_f U_f]')

err = abs( (I_f-U_f)/U_f );

if err>1e-9
    err
    error('Test:NetFlow','Compute net flow failed.');
    status = 0;
else
    disp('Compute net flow OK');
    status = 1;
end

return

% % Optional stuff for plotting result
% wmax = max(abs(W(:)));
% Xp=X;
% Yp=Y;
% Wp=W;

% [px py] = meshgrid(0:1);
% px=px(:);
% py=py(:);

% figure(1), clf
% for i=1:numel(px)
%     rx=px(i)*box(1);
%     ry=py(i)*box(2);
%     surf(Xp+rx, Yp+ry, Wp)
%     hold on
% end
% caxis([-wmax wmax])
% colorbar

% figure(2), clf
% for i=1:numel(px)
%     rx=px(i)*box(1);
%     ry=py(i)*box(2);
%     plot3(xvec(1,1)+rx, xvec(1,2)+ry, xvec(1,3),'or')
%     hold on
%     pcolor(Xp+rx,Yp+ry,Wp)
% end
% shading interp
% grid on
% zlim([0 box(3)])
% caxis([-wmax wmax])
% colorbar

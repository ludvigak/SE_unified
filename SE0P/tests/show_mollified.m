clear

box = [1 1 1];
R = sqrt(3);
M = 7*[1 1 1];
s = 4;

% Harmonic
solution = @(x) 1./abs(x);
exact_kernel = @(K,R) 4*pi*K.^-2;
kernel = @kernels.harmonic;

[k1 k2 k3] = k_vectors(s*M, s*box);
[K1 K2 K3] = ndgrid(k1,k2,k3);
K = sqrt(K1.^2 + K2.^2 + K3.^2);
GRhat = kernel(K, R);
GR = ifftn(ifftshift(GRhat));
% GR now has rubbish in the center
% truncate in real space by shifting, picking out center and shifting back
% usage of fftshift/ifftshift might not be right, but grid size is even anyway
n = M;
N = s*M;
GRc = ifftshift(GR);
GRt = GRc( floor((N(1)-n(1))/2) + (1:n(1)), ...
         floor((N(2)-n(2))/2) + (1:n(2)), ...
         floor((N(3)-n(3))/2) + (1:n(3)) );
GRts = fftshift(GRt);
GRthat = fftshift(fftn(GRts));


smid = (s*M(1))/2+1;
mid = (M(1)+1)/2;

kf = linspace(min(k1),max(k1),200);
[k1o k2o k3o] = k_vectors(M, box);


sfigure(1);clf
plot(kf, exact_kernel(kf), 'k--', 'DisplayName','$\hat{H}(k)$');
hold on
plot(kf, kernel(kf,R), '-', 'DisplayName','$\hat{H}^{\mathcal R}(k)$');
plot(k1o,GRthat(:,mid,mid), '.-', ...
     'DisplayName','$\hat{\tilde H}^{\mathcal R}(k)$')
xlabel('$k$','interpreter','latex')
axis([-10 10 0 10])
%grid on
h = legend('toggle');
set(h,'interpreter','latex')

sfigure(2);clf
h = box(1)/M(1);
x = h*((-M(1)/2):(M(1)/2-1) );
xs = h*((-s*M(1)/2):(s*M(1)/2-1) );
xf = [-2 -R-eps linspace(-R+eps,R-eps,100) R+eps 2];
volfunc = GRc/h^3;
Gtrunc = solution(xf); 
Gtrunc(abs(xf) > R) = 0;
plot(xf,Gtrunc, '--k', 'DisplayName','$H(r)$');
hold on
plot(xs,volfunc(:,smid,smid),'.-r', 'DisplayName','$\tilde H^{\mathcal R}(r)$');
axis([-1 1 0 20])
%grid on
xlabel('r')

% Plot Diagonal entries
xdia = [];
Gdia = [];
for i=1:numel(xs)
    xdia(i) = sqrt(3)*xs(i);
    Gdia(i) = volfunc(i,i,i);
end
xdia = sqrt(3)*xs;
Gdia = arrayfun(@(i) volfunc(i,i,i), 1:size(volfunc,1));
%plot(xdia,Gdia,'.-' , 'DisplayName','$\tilde{ H}^{\mathcal R}$')
h = legend('toggle');
set(h,'interpreter','latex')

input('ENTER to overwrite current plots, Ctrl-c to cancel');
write_tikz(1, ['../latex/fig/mollified_k'])
write_tikz(2, ['../latex/fig/mollified_r'])

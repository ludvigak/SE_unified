clear

filename_state = 'mfile_tests/refmat/test_error_kmax_ref_state.mat';
load(filename_state);

% % ## test integration of Bdir
% for k=[1 2 3 4]*10
%     for r=[1 3 5 7]
% I = zeros(3,3);
% for j=1:3
%     for l=1:3
%         for m=1:3
%             Iex(j,l,m) = Bdir_integral( j,l,m,k,r);
%             I(j,l,m) = ...
%                 integral2( @(theta, phi) Bdir(j,l,m,theta,phi,k,r).*sin(theta), ...
%                 0, pi, 0, 2*pi);
%         end
%     end
% end
% errinf = norm(Iex(:)-I(:),inf)
%     end
% end
% % ## works satisfactory

% % Visualize k integrand
% j=1
% l=1
% m=3
% k = linspace(1,30,10000);
% r = 2;
% xi = 3;
% I = kintegrand(j,l,m,k,r,xi);
% figure(1),clf
% plot(k,I)

xi = 7;
L = box(1);
kmax = 1:1:24;
kc = 2*kmax*pi/L;

Gsqint333 = zeros(size(kc));
rK =L/2;
Gsqint113 = zeros(size(kc));
for i=1:numel(kc)
    Gsqint333(i)=integral( @(r) (Gfunc(3,3,3,kc(i),r,xi)).^1.*r.^2, ...
                        0, rK, 'AbsTol',1e-3,'RelTol',1e-3 );
    Gsqint113(i)=integral( @(r) (Gfunc(1,1,3,kc(i),r,xi)).^1.*r.^2, ...
                        0, rK, 'AbsTol',1e-3,'RelTol',1e-3 );
end

%%
filename_state = 'mfile_tests/refmat/test_error_kmax_ref_state.mat';
load(filename_state);

S2 = zeros(3,3);
for l=1:3
    for m=1:3
        S2(l,m) = sum( ( f(:,l).*nvec(:,m) ).^2 );
    end
end
Q = sum(S2(:));
V = prod(box);

clf
erms = sqrt(Q)*(abs(Gsqint113)+abs(Gsqint333));

filename_errs = sprintf('mfile_tests/refmat/test_error_kmax_errs_xi=%.3f.mat',xi); 
load(filename_errs)
kmax_errs = 1:numel(errxi);

semilogy(kmax,erms, kmax, erms*errxi(1)/erms(1),'--'),hold all
semilogy(kmax_errs, errxi,'-*')

Bk = @(k,xi) exp(-k.^2/4/xi^2).*(8+2*k.^2/xi^2+k.^4/xi^4).*k/pi/V;

kmax2 = 1:30;
kc2 = 2*kmax2*pi/L;

semilogy(kmax2, FullErrFunc( kc2, L, xi ) * sqrt(Q), 'kd-')
11
x = linspace(1,30,200);
kx = 2*x*pi/L;
y = sqrt(Q)*exp(-kx.^2/4/xi^2).*sqrt(...
    321*sqrt(pi/2)*xi*erfc(kx/sqrt(2)/xi) + ...
    kx/xi^6.*( 257*xi^6 + 75*xi^4*kx.^2 + 11*xi^2*kx.^4 + kx.^6 ) ...
    )/10;

semilogy(x,y,'.')


ylim([1e-12 inf])
grid on
drawnow


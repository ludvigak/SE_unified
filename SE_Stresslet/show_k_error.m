clear

L = 1;
xi = 20;
N = 500;
idx = 1:N;
box = L*[1 1 1];
[x f n] = generate_state(N, box);
M_max = round(15*xi/pi);
M_ref = round(16*xi/pi);
opt.P = 32;
opt.box = box;
opt.M = M_ref*box;
ref = SE_Stresslet(idx, x, f, n, xi, opt);
ref_max = norm(ref(:), inf);

Mlist = 2:2:M_max;
err_inf = [];
err_rms = [];
est = [];
for i=1:numel(Mlist)
    this_opt = opt;
    this_opt.M = Mlist(i)*box;    
    SE_static  = SE_Stresslet_pre(x,xi,this_opt);
    uk = SE_Stresslet(idx, x, f, n, xi, this_opt, SE_static);
    err = uk - ref;
    err_rms(i) = sqrt(1/(3*N)*sum(err(:).^2));
end

% ESTIMATE
S2 = zeros(3,3);
for l=1:3
    for m=1:3
        S2(l,m) = sum( ( f(:,l).*n(:,m) ).^2 );
    end
end
Q = sum(S2(:));
V = prod(box);
L = min(box);
kmax = ( (Mlist-1)/2 );
kc = Mlist*pi;
kcxi = kc/xi;
est = sqrt(Q/V)*sqrt(L)*xi^1*exp(-kcxi.^2/4).*exp(.43*kcxi);



K = Mlist*pi;
clf
semilogy(K / xi, err_rms / ref_max, '.')    
hold on
semilogy(K / xi, est / ref_max,'-')
ylim([1e-16 1])
grid on
xlabel('$K/\xi$','interpreter','latex')
ylabel('Fourier space RMS error (rel.)')
legend('Measured','Estimate');
drawnow

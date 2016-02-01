clc; clear all
gamma_ = 0.577215664901532860606512090082402431042;

P = 10;
scale = 4;

z = [.26 .28];
z = [.1 .1];
z0 = z;

% mpole z with zj
center = [-.375 -.375];
a = compute_coeffs(P,z,scale,center);

% cluster around 0
N = 3;
rng(1);
zj = rand(N,2)-.5;
zj=[-.4 -.3 
    -.37 -.27
    -.42 -.33
    ];

q = rand(N+1,1);
q = q -mean(q);
qj = q(1:N);
mk = zeros(P+1,P+1);
for k1=0:P
    xk = (-zj(:,1)+center(1)).^k1;
    for k2=0:P
        mk(k1+1,k2+1) = sum(qj.*xk.*(-zj(:,2)+center(2)).^k2);
    end
end
for k1=0:P
    for k2=0:P
        S(k1+1,k2+1) = scale^(k1+k2);
    end
end

g(4) = sum(sum(mk.*a.*S));


% mpole zj with z
center = [.125 .125];
g(1:3) = 0;
for k1=0:P
    xk = (-z(1)+center(1))^k1;
    for k2=0:P
        mk(k1+1,k2+1) = q(end)*xk*(-z(2)+center(2))^k2;
    end
end

for k1=0:P
    for k2=0:P
        S(k1+1,k2+1) = scale^(k1+k2);
    end
end
for j=1:3
    a = compute_coeffs(P,zj(j,:),scale,center);
    g(j) = sum(sum(mk.*a.*S));
end

% direct
for j=1:3
    s = 0;
    for k=1:3
        if(j==k)
            continue;
        else
            x2y2 = (zj(j,1)-zj(k,1)).^2+(zj(j,2)-zj(k,2)).^2;
            s = s + q(k).*(expint(scale^2*x2y2)+log(scale^2*x2y2)+gamma_);
        end
    end
    g(j) = g(j) +s ;
end




% ak = a.*(abs(a)>1e-15);
% subplot(211)
% mesh(log(abs(ak)))
% axis([0 40 0 40 -15 0])
% subplot(212)
% mkk = mk.*(abs(mk)>1e-15);
% mesh(log(abs(mkk)))
% axis([0 40 0 40 -15 1])

% direct
zj = [zj; z];
for j=1:4
    s = 0;
    for k=1:4
        if(j==k)
            continue;
        else
            x2y2 = (zj(j,1)-zj(k,1)).^2+(zj(j,2)-zj(k,2)).^2;
            s =  s + q(k).*(expint(scale^2*x2y2)+log(scale^2*x2y2)+gamma_);
        end
    end
    f(j) = s;
end

err =sqrt(sum((f-g).^2)/sum(f.^2));
fprintf('%.16g\n',err);

return
mesh(X,Y,log(abs(a)))
xlabel('x')
ylabel('y')
hold on
%mesh(X,Y,log(abs(b)))

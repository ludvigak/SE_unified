
clear all, close all

box = [1 1 1];
M = 31*[1 1 1];
P = 19;

NN = [10 100 1000 2000];
c = 10;

opt.box=box;
opt.M=M;
opt.P=P;
opt.c=c;

for k = 1:length(NN)
    N = NN(k);
    x = rand(N,3);
    q = rand(N,1);

    % test FG gridding
    h1=SE_grid(x,q,box,M,P,c);
    h2=SE_fg_grid_mex(x,q,opt);

    max(abs(h1(:)-h2(:)))

    % test FG "integration"
    H = h1;
    
    phi_1 = zeros(N,1);
    for i = 1:N
        phi_1(i)=SE_int(x(i,:),H,box,M,P,c);
    end
    
    phi_2=SE_fg_int_mex(x,H,opt); 
    max(abs(phi_1-phi_2))
end


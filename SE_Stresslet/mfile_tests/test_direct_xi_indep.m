function status = test_direct_xi_indep()

format long e
rand('state',1)


N=10; 
TOL_r=1e-20;
NOL_r=12;
NOL_k=12; 

box = [1.1 1+pi/10 sqrt(2)/2];
[x f nvec] = generate_state(N,box);

xi = [1 2 5];

nxi=length(xi);
idx = 1;
ur=zeros(nxi,3);
uk=zeros(nxi,3);
utot=zeros(nxi,3);

tic
clf, leglist={};
for i = 1:nxi
    uk(i,:)=stresslet_direct_fd  (idx, x, f, nvec, xi(i),  box, NOL_k)';
    [rtemp conv]=stresslet_direct_real(idx, x, f, nvec, xi(i),  box, NOL_r,TOL_r);
    ur(i,:) = rtemp';
    utot(i,:)=ur(i,:)+uk(i,:);
    %%No self interaction term for stresslet. 
    
    disp(['idx= ' num2str(idx) ', xi=' num2str(xi(i)) '.']); 
    disp('ur= ');
    disp(ur(i,:))
    disp('uk ');
    disp(uk(i,:))
    disp('utot');
    disp(utot(i,:))
    
    semilogy(conv,'.-'), hold all
    leglist{end+1}=['\xi=' num2str(xi(i))];
end;

legend(leglist), ylabel('Real space layer contrib.'), xlabel('Layer'), grid on
toc

utot
xi
box

res=abs(diff(utot))

if max(res(:))<1e-10
    status = 1;
    fprintf('\n********** EWALD XI INDEPENDENCE: OK **********\n\n')
else
    status = 0;
    warning('EWALD XI INDEPENDENCE: FAILED')
end


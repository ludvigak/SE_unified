%% figure 8 example 2 in paper 
%% total runtime in potential 
%% with 1e-6 accuracy
%% with icc we get better than linear!!!

clear
rng(1);

Llist = 6:14;
Nlist = 100*Llist.^3;
Mlist = [46 54 62 64 80 84 92 100 108];
NLlist = [4 5 6 6 6 7 8 9 10];

time = [];
for k=1:numel(Llist)
    N = Nlist(k);
    L = Llist(k);
    box = [L L L];
    [x, f] = vector_system(N, box);

    opt.M = Mlist(k);
    opt.xi = 3.5;
    opt.P = 16;
    opt.rc = .9;
    opt.box = box;
    opt.layers = 5;
    opt.sl = 3;
    opt.nl = NLlist(k);
    opt.s0 = 2.5;

    
    %% Ewald
    disp('One periodic Ewald...')
    [uf tfourier]= se1p_fourier_space(x, f, opt);
    tr=tic;
    ur = SE1P_rsrc_cell_mex(x', f, opt.rc, opt.xi, opt.layers, box);
    treal=toc(tr);
    ue = uf+ur;
    time(end+1) = tfourier.total+treal;
end
plot(Nlist,time,'b.')
hold on
plot(Nlist,Nlist./Nlist(1)*time(1),'r--')

%% Direct
% idx = 1:N;
% disp('Direct sum...')
% opt.P=24;
% u1 = se1p_fourier_space_direct(idx,x, f, opt, false);
% u2 = se1p_real_space(idx,x, f, opt,true);
% u3 = se1p_k0_direct(idx,x,f,opt,true);
% u  = u1+u2+u3;

% ue = ue(idx,:);
% rms_error = rms(u(:)-ue(:)) / rms(u(:))

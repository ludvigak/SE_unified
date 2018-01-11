function shape_parameter 
    global N
    
restoredefaultpath
addpath('..')

init

rng(1)
box = [1 1 1];   % domain
N = 1e2;          % number of charged particles
M0 = 28;

opt.M = M0*box;
opt.P = 20;
opt.xi = pi*M0 / 12; opt.xi = 6.3;
opt.box = box;
opt.layers = (opt.M(1)-1)/2;
opt.beta = 2.5;

% charge-neutral system
[x q] = vector_system(N,box);

Plist = 4:2:24;
Plist0 = 2:2:18;
colk={'ro-','bs-','g^-','k'};
colg={'ko-','ks-','kd-','k>-'};

%% --------------------------- Iterate Kaiser ------------------------------
%% 3d-periodic
opt_3p = opt;
disp('direct 3p')
ref3p = compute_direct_3p(x,q,opt_3p);
sfigure(1);
clf
hold on
i=1;
for b = [15 25 30 1]
    disp('Kaiser 3p')
    for k=1:numel(Plist0)
        opt_3p.P = Plist0(k);
        if(b~=1)
            opt_3p.beta = b/Plist0(k);
        else
            opt_3p.beta = 2.5;
        end
        [~, kaiser(k)] = call_kaiser_3p(x, q, opt_3p, ref3p);
    end
    plot(Plist0,kaiser,colk{i},'DisplayName',['\beta=', num2str(b)])
    i = i + 1;
end

plot(Plist0,10*exp(-2.5*Plist0),'k--')    
set(gca,'yscale','log')
xlabel('P')
ylabel('e_{rms} (rel.)')
h = legend('toggle');
set(h, 'Location','SouthWest')
ylim([1e-15 1e0])
xlim([4 14])


%% ------------------------------ Iterate gaussian --------------------------
%% 3d-periodic
disp('gaussian 3p')
opt_3p = opt;
addpath('../SE')
sfigure(2)
clf
hold on
i=1;
for m=[5 6 7 1]
    for k=1:numel(Plist)
        opt_3p.P = Plist(k);
        if(m~=1)
            opt_3p.m = m;
        else
            opt_3p.m=.95*sqrt(pi*Plist(k));
        end        
        [~, gaussian(k)] = call_gaussian_3p(x, q, opt_3p, ref3p);
    end
    plot(Plist,gaussian,colk{i},'DisplayName',['m=', num2str(opt_3p.m)])
    i = i + 1;
end

%% -------------------------------- PLOT ------------------------------

set(gca,'yscale','log')
xlabel('P')
ylabel('e_{rms} (rel.)')
h = legend('toggle');
set(h, 'Location','SouthWest')
ylim([1e-16 1e0])

end

%%% ---------------------------- COMPUTE DIRECT ------------------------------
%% 3d-periodic direct
function ref3p = compute_direct_3p(x,q,opt)
    global N
% idx = 1:10;
%ref3p = SE3P_direct_fd_mex(idx,x,q,opt);
    ref3p = se3p_fourier_space_kaiser(1:N,x,q,opt);
end

%%% ---------------------------- CALL KAISER ------------------------------
function [u kaiser_3p] = call_kaiser_3p(x, q, opt, ref3p)
    global NOPRE N

    [u time] = se3p_fourier_space_kaiser(1:N,x,q,opt);
    kaiser_3p = rms(u-ref3p)/rms(ref3p);
end

%%% ---------------------------- CALL GAUSSIAN ------------------------------
function [u gaussian_3p] = ...
        call_gaussian_3p(x, q, opt, ref3p)
    
    global NOPRE
    u   = se3p_fourier_space(x,q,opt);
    gaussian_3p   = rms(u-ref3p)/rms(ref3p);
end
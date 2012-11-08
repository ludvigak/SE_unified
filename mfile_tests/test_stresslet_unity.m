function status = test_stresslet_unity()

% Test integrating unity constant stresslet distribution over closed body

clear


inside = 0;
onsurf = 0;
    
L       = 1;
nbox    = 1;
xi      = 5;
R       = 0.25;

% non zero component
jj = 3;

SE_opt.box = [L L L];
SE_opt.M = 16*[1 1 1];
SE_opt.P = 16;  
SE_opt.xi = xi;

Nlist = [10 20 40 80 160];
err = [];

for iN = 1:numel(Nlist)
    N = Nlist(iN);

    bdys = 1;

    % Nystrom disc., trapezoidal rule
    % Punctured trapezoidal on surface
    
    theta = linspace(0,2*pi,N+1);
    wtheta = diff(theta);
    theta = theta(1:end-1);

    phi = linspace(0,pi,N+1)';
    wphi = diff(phi);
    phi = phi(1:end-1)+wphi/2;


    [THETA PHI] = meshgrid(theta,phi);

    X = R*cos(THETA).*sin(PHI);
    Y = R*sin(THETA).*sin(PHI);
    Z = R*cos(PHI);

    NX = X/R;
    NY = Y/R;
    NZ = Z/R;

    X0 = [0.5 0.5 0.5];

    X = X+X0(1);
    Y = Y+X0(2);
    Z = Z+X0(3);

    Q = (R^2*sin(phi).*wphi) * wtheta;

    xvec = [X(:) Y(:) Z(:)];
    nvec = [NX(:) NY(:) NZ(:)];
    fvec = zeros(size(xvec));
    fvec(:,jj)=1;
    qweight = Q(:);
    nvec = bsxfun(@times, nvec, qweight);
    
    
    % measure points outside, on, inside
    Ifree = [0 4*pi 8*pi];
    xs = [ [0 0 0]+pi/100; X0];
    ns = [ 0 0 0; 0 0 0];
    fs = [ 0 0 0; 0 0 0];
    
    fvec = [fs; fvec];
    xvec = [xs; xvec];
    nvec = [ns; nvec];
    
    idx = [1 N/2+2 2];

    SE_static  = SE_Stresslet_pre(xvec,xi,SE_opt);
    
    ewres_r = stresslet_direct_real_fast( idx, xvec, fvec, nvec, xi, SE_opt.box, nbox, 1e100);
    ewres_f = SE_Stresslet(idx,xvec,fvec,nvec,xi,SE_opt,SE_static);
    ewres_z = stresslet_direct_fd_zero( idx, xvec, fvec, nvec, SE_opt.box);

    ewres = ewres_f + ewres_r + ewres_z;
    err = [err; ewres(:,jj)'-Ifree]; % +8*pi*conc
    
end
%%
format long
err=abs(err)
k = [0 0 0];
for i=1:3
    c = polyfit(-log(Nlist),log(err(:,i)'),1);
    k(i)=c(1);
end
format short
k
err_k = abs(k - [2 1 2])
if any(err_k>0.05)
    warning('Test:StressletUnity','Integration of stresslet unity failed.');
    status = 0;
else
    disp('Integration of stresslet unity OK.');
    status = 1;
end
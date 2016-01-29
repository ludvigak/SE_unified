function main(P, N)
    addpath('bin');
    addpath('../SE_Stresslet')

    
    rng(1)
    
    for P=[P]    
        clear output;
        box = [1 1 1];       
        x = bsxfun(@times, rand(N, 3), box);
        f = rand(N, 1);   
        params.box = box;
        params.M = box*100;
        params.P = P;
        params.xi = 10;
        opt = setup(params);
        [zx zy zz idx] = SE_fgg_expand_all_mex(x,opt);
        zs = SE_fgg_base_gaussian_mex(opt);
        
        formatter = @(n,e,t) struct('name', n,'error', e,'time', t, 'rel_time', 0);
        % Vanilla / ref
        F_ref = SE_fg_grid_mex(x,f,opt);
        t_ref = timeit(@() SE_fg_grid_mex(x,f,opt) );       
        output(1) = formatter('MEX FGG / REF',NaN,t_ref);
        relerr = @(X) norm(F_ref(:)-X(:), inf) / norm(F_ref(:), inf);
        
        % MEX SSE FGG
        F_fgg = SE_fg_grid_split_mex(x,f,opt,zs,zx,zy,zz,idx);
        t_fgg = timeit(@() SE_fg_grid_split_mex(x,f,opt,zs,zx,zy,zz,idx) );    
        output(end+1) = formatter('MEX FGG AVX2', relerr(F_fgg), t_fgg);

        % THRD
        F_fggp = SE_fg_grid_thrd_mex(x,f,opt);
        t_fggp = timeit(@() SE_fg_grid_thrd_mex(x,f,opt) );    
        output(end+1) = formatter('MEX FGG THRD', relerr(F_fggp), t_fggp);

        % THRD SPLIT
        F_fgg = SE_fg_grid_split_thrd_mex(x,f,opt,zs,zx,zy,zz,idx);
        t_fgg = timeit(@() SE_fg_grid_split_thrd_mex(x,f,opt,zs,zx,zy,zz,idx) );    
        output(end+1) = formatter('MEX FGG AVX2 THRD', relerr(F_fgg), t_fgg);        
        
        for i=1:numel(output)
            output(i).rel_time = output(i).time / output(1).time;
        end
        
        % vanilla
        % s = tic();
        % F_van = vanilla_grid(x, f, params.box, params.M, params.P, opt.c);
        % t_van = toc(s);    
        % output(end+1) = formatter('Vanilla', relerr(F_van), t_van);    
        
        % Inverted vanilla
        % s = tic();
        % F_inv = inv_grid(x, f, params.box, params.M, params.P, opt.c);
        % t_inv = toc(s);    
        % output(end+1) = formatter('Inverted', relerr(F_inv), t_inv);    

        
        
        % Output
        disp(' ');
        disp('Run parameters:')
        fprintf('N = %d\nP = %d \n', N, P);
        fprintf('M = [%s]\n', num2str(params.M));
        disp(struct2table(output))
        if max([output.error]) > 1e-10
            error('FAILURE!')
        end
    end
end

function H=vanilla_grid(xvec, fvec, box, M, P, c)
    H = zeros(M);
    N = size(xvec,1);
    h = box(1)/M(1);
    for i=1:N  
        x = xvec(i,:);
        f = fvec(i,:);
       idx = round(x/h)+1;
       idx_low = idx-(P-1)/2;
       idx_upp = idx+(P-1)/2;       
       for supp_x = idx_low(1):idx_upp(1);
           for supp_y = idx_low(2):idx_upp(2);
               for supp_z = idx_low(3):idx_upp(3);
                   
                   % grid points (assume interval start at 0)
                   PX = h*(supp_x-1);
                   PY = h*(supp_y-1);
                   PZ = h*(supp_z-1);

                   % distances
                   D2 = (PX-x(1)).^2 + (PY-x(2)).^2 + (PZ-x(3)).^2;

                   % assignment indices (wrapped to produce periodicity)
                   ix = mod(supp_x-1,M(1))+1;
                   iy = mod(supp_y-1,M(2))+1;
                   iz = mod(supp_z-1,M(3))+1;

                   % gaussian
                   Q = (c/pi)^1.5*exp(-c*D2);
                   H(ix, iy, iz) = H(ix, iy, iz) + Q*f;
               end
           end
       end        
    end
end

function H=inv_grid(xvec, fvec, box, M, P, c)
    H = zeros(M);
    N = size(xvec,1);
    h = box(1)/M(1);
    rc2 = (h*(P-1)/2)^2;
    parfor iz=1:M(3)
        Hz = zeros(M(1),M(2));
        for iy=1:M(2)                       
            for ix=1:M(1)          
                val = 0;
                idx_grid = [ix iy iz];
                pos = h*(idx_grid-1);
                for i=1:N
                    x = xvec(i,:);
                    f = fvec(i);
                    idx_x = round(x/h)+1;
                    
                    %p = closest_image(idx_grid, idx_x, M);                    
                    % Because MATLAB is crap:
                    p = [0 0 0];
                    for ip=1:3
                        gd = abs(idx_grid(ip)+M(ip)*(-1:1)-idx_x(ip));
                        if gd(1) < gd(2)
                            if gd(1) < gd(3)
                                p(ip) = gd(1);
                            else
                                p(ip) = gd(3);
                            end
                        else
                            if gd(2) < gd(3)
                                p(ip) = gd(2);
                            else
                                p(ip) = gd(3);
                            end
                        end
                    end
                                        
                    % idx_grid_p = idx_grid + p.*M;
                    % d_grid = idx_x - idx_grid_p;
                    % if max(abs(d_grid)) > (P-1)/2
                    %     continue
                    % end
                    d = x - (pos + p.*box);
                    d2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3);
                    if d2 <= rc2
                        val = val + exp(-c*d2)*f;
                    end
                end
                Hz(ix, iy) = (c/pi)^1.5*val;
            end
        end
        H(:,:,iz) = Hz; 
    end
end

function shift_t = closest_image(idx_t, idx_s, M)
    shift_t = [0 0 0];
    for i=1:3
        [~, min_comp] = min(abs(idx_t(i)+M(i)*(-1:1)-idx_s(i)));
        shift_t(i) = min_comp-2;
    end
end


function opt = setup(params)   
% parameters and constants
    xi = params.xi;
    opt = parse_params(params);
    [w m] = unpack_params(opt);
    eta = (2*w*xi/m)^2;
    opt.c = 2*xi^2/eta;
end


function [w m M P] = unpack_params(opt)
    w = opt.w;
    m = opt.m;
    M = opt.M;
    P = opt.P;
end

function [H info] = SE_grid(x, q, box, M, P, c)
    % grid size
    h = box(1)/M(1);

    % input checking
    assert(abs( h - box(2)/M(2)) < eps &&  abs( h - box(3)/M(3)) < eps);
    assert(P <= min(M))
    assert(mod(P-1,2) < eps)

    H = zeros(M);
    N = size(x,1);

    for i=1:N
        [Q ax ay az] = SE_gaussian(x(i,:),h,M,P,c);
        
        % put into grid fcn
        H(ax,ay,az) = H(ax,ay,az) + q(i)*Q;
    end

    info.gaussian_mass_resid = abs(1-h^3*sum(Q(:)));
    info.arithmetic_ratio = P^3/prod(M);
    info.h = h;
end

function [Q ix iy iz] = SE_gaussian(x,h,M,P,c)
   % support indices
    idx = round(x/h)+1;
    idx_low = idx-(P-1)/2;
    idx_upp = idx+(P-1)/2;

    supp_x = idx_low(1):idx_upp(1);
    supp_y = idx_low(2):idx_upp(2);
    supp_z = idx_low(3):idx_upp(3);

    [SX SY SZ] = ndgrid(supp_x,supp_y,supp_z);

    % grid points (assume interval start at 0)
    PX = h*(SX-1);
    PY = h*(SY-1);
    PZ = h*(SZ-1);

    % distances
    D2 = (PX-x(1)).^2 + (PY-x(2)).^2 + (PZ-x(3)).^2;

    % assignment indices (wrapped to produce periodicity)
    ix = mod(supp_x-1,M(1))+1;
    iy = mod(supp_y-1,M(2))+1;
    iz = mod(supp_z-1,M(3))+1;

    % gaussian
    Q = (c/pi)^1.5*exp(-c*D2);
end

function fgg_timing(P, N)
    rng('default')
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
        output(1) = formatter('MEX FGG (reference)',0,t_ref);
        relerr = @(X) norm(F_ref(:)-X(:), inf) / norm(F_ref(:), inf);
        
        % MEX SSE FGG
        F_fgg = SE_fg_grid_split_mex(x,f,opt,zs,zx,zy,zz,idx);
        t_fgg = timeit(@() SE_fg_grid_split_mex(x,f,opt,zs,zx,zy,zz,idx) );    
        output(end+1) = formatter('MEX FGG SIMD', relerr(F_fgg), t_fgg);

        % THRD
        F_fggp = SE_fg_grid_thrd_mex(x,f,opt);
        t_fggp = timeit(@() SE_fg_grid_thrd_mex(x,f,opt) );    
        output(end+1) = formatter('MEX FGG THRD', relerr(F_fggp), t_fggp);

        % THRD SPLIT
        F_fgg = SE_fg_grid_split_thrd_mex(x,f,opt,zs,zx,zy,zz,idx);
        t_fgg = timeit(@() SE_fg_grid_split_thrd_mex(x,f,opt,zs,zx,zy,zz,idx) );    
        output(end+1) = formatter('MEX FGG SIMD THRD', relerr(F_fgg), t_fgg);        
        
        if mod(P,2) && N <= 1000
            % MATLAB
            s = tic();
            F_van = SE_grid(x, f, params.box, params.M, params.P, opt.c);
            t_van = toc(s);    
            output(end+1) = formatter('MATLAB FGG (gold standard)', relerr(F_van), t_van);
        end

        for i=1:numel(output)
            output(i).rel_time = output(i).time / output(1).time;
        end        
        
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

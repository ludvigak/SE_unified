clear all, close all,

test = 'a'

switch test
    
    case 'a'
        L = 3; box = [L L L];
        N = 13500;
        
        Nc = 1000
        
        (4*pi*N/(3*Nc*L^3))^(1/3)
        
        S = 9;
        
        xi = 10;
        
        [x q] = SE_state(N, box);
        %u = fast_real_sum(1:N, x, f, xi, @hasimoto_op_real, L, S);
        %tic
        [u wt] = fast_rs_sum(x, q, xi, L, S, true);
        disp(wt)
        %toc
        
    case 'b'
        
        L = 1; box = [L L L];
        N = 1000;
        S = 3;
        xi = 15;
        [x f] = generate_state(N, box, 1);
        x(1,:) = [1/2 1/2 1/2]; 
        
        format long 
        u1 = stokes_ewald_direct_real_mex(1, x', f', xi, box, 0)'
        u1 = stokes_ewald_direct_real_mex(1, x', f', xi, box, 1)'
        
        u = fast_rs_sum_II(x, f, xi, L, S, false);
        u(1,:)
        
end
clear

xilist = [6 8 10 14];

generate_ref_state = false;
generate_ref_solutions = generate_ref_state;

filename_state = 'mfile_tests/test_error_rc_ref_state.mat';
if generate_ref_state
    disp('Generating reference state...')
    N = 3000;
    box = [2 2 2];
    [x f nvec] = generate_state(N,box);
    idx = 1:N;
    save(filename_state,...
          'N','x','f','nvec','box','idx');
    disp('done.')
    fprintf('Saved %s\n',filename_state);
else
    load(filename_state);
end

if generate_ref_solutions
    for xi=xilist
        filename_sol = sprintf('mfile_tests/test_error_rc_ref_sol_xi=%.3f.mat',xi);    
        NOL_ref = 1;
        rc_ref = NOL_ref*max(box)*10;
        UREF = stresslet_direct_real_fast(idx, x, f, nvec, xi,  box, NOL_ref, rc_ref);
        save(filename_sol,...
              'N','x','f','nvec','box','idx','UREF');
        fprintf('Saved %s\n',filename_sol);
    end
end

figure(1),clf
figure(2),clf

kvot = [];

for xi=xilist
    filename_sol = sprintf('mfile_tests/test_error_rc_ref_sol_xi=%.3f.mat',xi);  
    load(filename_sol)
    NOL = 1;
    rc = 0.1:0.02:1;

    for i=1:numel(rc)
        if rc(i)>min(box)/2
            u = stresslet_direct_real_fast(idx, x, f, nvec, xi,  box, NOL, rc(i) );
        else
            u = stresslet_real_rc( x, f, nvec, xi, box, rc(i));
        end
        e = (u - UREF).^2; 
        err(i,:) = sqrt( sum(e(:))/(3*N) );
    end

    S2 = zeros(3,3);
    Qt = 0;
    for l=1:3
        for m=1:3
            S2(l,m) = sum( ( f(:,l).*nvec(:,m) ).^2 );
        end
    end
   
    Q = sum(S2(:))/9;
    V = prod(box);
        
    % Power of Mathematica!
    est = sqrt(Q/V * exp(-2*xi^2*rc.^2).*(...
            109*xi^2*rc + 1588/3*rc.^3*xi^4 -464*rc.^5*xi^6 + 448/3*rc.^7*xi^8 ...
          ) );
      
    % Short form
%     est = 8*sqrt(7*Q/V/3*rc.^7)*xi^4.*exp(-xi^2*rc.^2);
          
    figure(1)
    semilogy(rc, err,'.-', rc, est, '--r')  
    hold on
    ylim([1e-14 inf]);
    
    figure(2)
    thiskvot = est'./err;
    plot_idx = thiskvot>0.9;
    plot(rc(plot_idx), thiskvot(plot_idx),'.-')
    hold all
    
    kvot = [kvot; thiskvot];
end

plot(rc, rc.^0, '--k')

%%
figure(1)
title(['Truncation error of real space part (RMS), \xi = '...
        sprintf('%g,',xilist) ]);
xlabel('r_c')
ylabel('E_{RMS}')
ylim([1e-14 inf]);
grid on

figure(2)
leglist={};
for xi=xilist
    leglist{end+1} = sprintf('\\xi=%g',xi);
end
legend(leglist)
ylim([0.9 1.1])
ylabel('Estimate / Measured (RMS)')
xlabel('r_c')
title('Accuracy of real space truncation error estimate')
clear

if matlabpool('size')==0
    matlabpool('open')
end

test_direct_xi_indep
test_direct_real_fast
test_fast_summation
test_accuracy
test_sse

disp(' ')
disp(' ')
disp('**********************************')
disp('********** ALL TESTS OK **********')
disp('**********************************')
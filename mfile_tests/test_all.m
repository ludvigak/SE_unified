clear
close all

disp(' ')
disp(' ')
disp('***************************************')
disp('********** RUNNING ALL TESTS **********')
disp('***************************************')

if matlabpool('size')==0
    matlabpool('open')
end

test_direct_xi_indep
test_direct_real_fast
test_real_rc
test_rs_sing_sub
test_fast_summation
test_accuracy
test_sse

disp(' ')
disp(' ')
disp('**********************************')
disp('********** ALL TESTS OK **********')
disp('**********************************')
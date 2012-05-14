clear

if matlabpool('size')==0
    matlabpool('open')
end

test_direct_xi_indep
test_fast_summation
test_accuracy
test_sse
function test_all
clc
disp('Running all tests...')

if matlabpool('size')==0
    matlabpool('open')
end

d = dir('mfile_tests/test_*.m');
tests = {};
results = [];
for i=1:numel(d);
    fname = d(i).name(1:end-2);
    if strcmp(fname,'test_all')
        continue
    end
    tests{end+1} = fname;
    disp(' ');
    disp(['================================== '...
          fname ...
          ' ================================== ']);
    disp(' ');
    results(end+1) = feval(fname);
end

%%
disp(' ');
disp(['=================================='...
      '================' ...
      '================================== ']);
disp(' ')
disp('All tests completed:')
disp('--------------------')
N = numel(tests);
padder = '........................................';
for i=1:N
    fprintf('(%d/%d) %s%s', i, N, tests{i}, padder(1:(numel(padder)-numel(tests{i}))));
    if results(i)
        fprintf('ok\n')
    else
        fprintf('FAILED\n')
    end
end
disp('--------------------')
if all(results)
    disp('ALL TESTS OK')
else
    fprintf('%d/%d TESTS FAILED\n',numel(find(i>0)),N)
end
disp(' ')
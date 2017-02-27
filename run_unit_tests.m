clear

module_folders = {'SE_Rotlet', ...
                  'SE_Stresslet', ...
                  'SE_Stokes', ...
                  'SE_Stokes_direct'
                 };

% Run all unit tests that we have
results = cell(size(module_folders));
num_tests = numel(module_folders);
% Enter folders and run init + tests
for i=1:num_tests
    folder = module_folders{i};
    run([folder '/init'])
    results{i} = runtests([folder '/tests']);    
end
% Print results
all_passed = true;
for i=1:num_tests
    folder = module_folders{i};
    fprintf('======================================\n');
    fprintf('Results from %s:\n\n', folder);
    disp(table(results{i}))
    failed = find([results{i}.Passed] == false);
    if any(failed)
        all_passed = false;
    end
end
% List failed tests
if all_passed
    disp('ALL TESTS PASSED')
else
    disp('===================================')
    fprintf('= FAILED TESTS:\n=\n')
    for i=1:num_tests
        folder = module_folders{i};
        failed = find([results{i}.Passed] == false);
        if any(failed)
            for f=failed
                fprintf('= %s | %s\n', folder, results{i}(f).Name);
            end
        end
    end    
    disp('===================================')
    error('There were failed tests.')
end


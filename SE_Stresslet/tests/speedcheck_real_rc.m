near_n = 100;

M = 2:15;
rc = 1./M;
N = round(3*near_n./(4*pi*rc.^3));

% N = round( linspace(1000,300000,15) );
% N = [1000:1000:40000]
% rc = ( 3*near_n./(4*pi*N) ).^(1/3);

xi = 1;
box = [1 1 1];


figure(1), clf
for threads=1:4
    maxNumCompThreads(threads)

    [wall_time matrix_size] = deal(zeros(size(N)));
    for iN = 1:length(N)
        fprintf('----------------\nN=%d\n',N(iN));
        [x f nvec] = generate_state(N(iN),box);
        a = tic;
%         [res AMAT] = stresslet_real_rc( x, f, nvec, xi, box, rc(iN));
        res = stresslet_real_rc( x, f, nvec, xi, box, rc(iN));
        wall_time(iN) = toc(a);
%         matrix_size(iN) = getfield(whos('AMAT'), 'bytes');
%         fprintf('near_n_avg=%.1f\n',round(nnz(AMAT{1,1})/N(iN)))
        clear AMAT res;
    end

    %%

    plot(N, wall_time,'.-', N, N*wall_time(end)/N(end), 'k--'), hold all
    xlabel('N')
    ylabel('t [s]');
    legend('run time','O(N)','Location','NorthWest')
end

% figure(2), clf
% matrix_size_MB = matrix_size/1024^2;
% plot(N, matrix_size_MB,'.-', N, N*matrix_size_MB(end)/N(end), 'k--')
% xlabel('N')
% ylabel('MiB')
% legend('matrix size','O(N)');
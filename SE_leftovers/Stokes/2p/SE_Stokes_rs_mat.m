function [A stats] = SE_Stokes_rs_mat(x,rc,xi,varargin)
% Stokes Ewald real-space matrix assembler
% Dag Lindbo, dag@kth.se, Aug. 2011.
    
    N = size(x,1);
    verb = true;
    
    if rc<0
        % create full matrix for all-to-all         
        assert(N<=5000,'Problem probably too big')

        % initialize cell array to hold matirx blocks
        for i=1:3
            for j=1:3
                G{i,j} = zeros(N,N);
            end
        end
        
        % for each x_0
        for m=1:N
            
            % current point
            x0 = x(m,:);
            
            % distance to other points (incl. self)
            xh = x-repmat(x0,[N 1]);
            r2 = sum(xh.^2 ,2);
            r = sqrt(r2);
            
            % evaluate rs operator
            b = xi^2*r2;        
            g = 2*(xi*exp(-b)./(sqrt(pi)*r2) + erfc(xi*r)./(2*r.*r2));
            g(m) = 0;
            f = r2.*g - 4*xi*exp(-b)/sqrt(pi);
            f(m) = 0;        
            for i=1:3
                for j=1:3
                    G{i,j}(m,:) = f*(i==j) + g.*xh(:,i).*xh(:,j);
                end
            end
            
        end
        
        % concatenate to full operator
        A = cell2mat(G);
        
    else
        % assemble sparse matrix, which is equivalent to neighbour
        % list RS sum

        rowvec = [];
        colvec = [];
        Avec = [];
        %prog = progress_bar();
        
        if nargin == 3
            % no help
            cprintf(verb,'[MAKE RS MAT] (sparse) N=%d\t\n',N)
            
            for m=1:N
                %if mod(m,500)==0
                %    prog.update(m/N);
                %end
                
                % current point
                x0 = x(m,:);
                
                % distance to other points (incl. self)
                xh = [x(:,1)-x0(1) x(:,2)-x0(2) x(:,3)-x0(3)];
                r2 = sum(xh.^2 ,2);

                % find neighbours, and exclude self
                nbh_idx = find(r2<=rc^2 & r2>0);

                if(~isempty(nbh_idx))
                    
                    % remove non-neighbours
                    xh = xh(nbh_idx,:);
                    r2 = r2(nbh_idx);
                    
                    assert(all(r2>0)) % self should have been removed
                    
                    % evaluate operator (radial parts)
                    r = sqrt(r2);
                    b = xi^2*r2;        
                    g = 2*(xi*exp(-b)./(sqrt(pi)*r2) + erfc(xi*r)./(2*r.*r2));
                    f = r2.*g - 4*xi*exp(-b)/sqrt(pi);
                    
                    % evaluate operator and add to triplet vectors
                    row = m*ones(size(nbh_idx));
                    for i=1:3
                        im = row + N*(i-1);
                        for j=1:3
                            jm = nbh_idx + N*(j-1); 
                            a = f*(i==j) + g.*xh(:,i).*xh(:,j);
                            
                            rowvec = [rowvec; im];
                            colvec = [colvec; jm];
                            Avec = [Avec;a];
                        end
                    end
                end
            end
            
        else
            % use additional information about which dofs are close
            
            cprintf(verb,'[MAKE RS MAT] (sparse, nbh list) N=%d\t\n',N)
            prog_msg = ' Computing matrix elements on row %7d (%d)\r';
            nbh_blocks = varargin{1};
            block_sz = varargin{2};
            nbh_shift = varargin{3};
            e=1:block_sz;            
            assert(block_sz*length(nbh_blocks) == N)
            
            nnz_max = get_nnz_est(block_sz, nbh_blocks, verb);
            rowvec = zeros(nnz_max,1);
            colvec = zeros(nnz_max,1);
            Avec = zeros(nnz_max,1);
            pos = 0;
            
            tic()
            for m=1:N
                if mod(m,100)==0
                    cprintf(verb,prog_msg,m,N);
                end
                
                % current point
                x0 = x(m,:);
                block_idx = idivide(int32(m-1),int32(block_sz))+1;

                % build list of neighbours to search
                [nbh_x nbh_idx] = ...
                    get_nbh_pts(m, x, nbh_blocks, block_sz, nbh_shift);
                
                % remove self
                self = nbh_idx==m;
                nbh_idx(self) = [];
                nbh_x(  self,:) = [];

                if(~isempty(nbh_idx))
                    % distance to other points
                    xh = [nbh_x(:,1)-x0(1) nbh_x(:,2)-x0(2) nbh_x(:,3)-x0(3)];
                    r2 = sum(xh.^2 ,2);

                    % find neighbours, and exclude self
                    inside_rc = r2<=rc^2;
                    
                    % remove non-neighbours
                    nbh_idx = nbh_idx(inside_rc);
                    xh = xh(inside_rc,:);
                    r2 = r2(inside_rc);
                    
                    assert(all(r2>0)) % self should have been removed
                    
                    % evaluate operator (radial parts)
                    r = sqrt(r2);
                    b = xi^2*r2;        
                    g = 2*(xi*exp(-b)./(sqrt(pi)*r2) + erfc(xi*r)./(2*r.*r2));
                    f = r2.*g - 4*xi*exp(-b)/sqrt(pi);
                    
                    % evaluate operator and add to triplet vectors
                    row = m*ones(size(nbh_idx));
                    nz=length(row);
                    for i=1:3
                        % row indices
                        im = row + N*(i-1);
                        for j=1:3                            
                            % column indices
                            jm = nbh_idx + N*(j-1); 
                            
                            % matrix elements
                            a = f*(i==j) + g.*xh(:,i).*xh(:,j);

                            % insert into vectos
                            rowvec( (pos+1):(pos+nz)) = im;
                            colvec( (pos+1):(pos+nz)) = jm;
                            Avec(   (pos+1):(pos+nz)) = a;
                            pos = pos+nz;
                        end
                    end
                end
            end
        end
        stats.wtime(1) = toc();
        cprintf(verb,prog_msg,m,N);
        cprintf(verb,'\n')
        cprintf(verb,' Assembling.... ');
        tic
        A = sparse(rowvec(1:pos),colvec(1:pos),Avec(1:pos),3*N,3*N);
        stats.wtime(2)=toc();
        cprintf(verb,' Done. NNZ = %.2e, NNZ/row/comp = %.1f, fill = %.2f\n',...
                nnz(A), nnz(A)/(9*N),nnz(A)/prod(size(A)));
    end
    
function [p idx] = get_nbh_pts(m,x,nbh_blocks,block_sz,nbh_shift)

    e=1:block_sz;
    block_idx = idivide(int32(m-1),int32(block_sz))+1;
    
    % indices
    m_block_idx=nbh_blocks{block_idx};
    idx = zeros(block_sz*length(m_block_idx),1);
    for k=1:length(m_block_idx)
        idx(((k-1)*block_sz:(k*block_sz-1))+1) = ...
            e+(m_block_idx(k)-1)*block_sz;
    end
    
    % shifts to account for periodicity
    e = ones(block_sz,1);
    s = [];
    for k=1:3
        q = nbh_shift{block_idx}(:,k)';
        qq = q(e,:);
        s = [s qq(:)];
    end

    % points
    p = x(idx,:) + s;
    
function n = get_nnz_est(block_sz, nbh_blocks,verb)
    
    m = length(nbh_blocks);
    n = 0;
    for j=1:m
        nn = length(nbh_blocks{j});
        n = n+nn;
    end
    n = n*block_sz^2*9;
    cprintf(verb,' NNZ upper bound: %.2e \n',n)
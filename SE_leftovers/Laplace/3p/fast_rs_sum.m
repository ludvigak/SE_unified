function [phi wtime] = fast_rs_sum(x, q, xi, L, S, verlet)

plotting = false;
checking = false;
verbose = true;

N = size(x,1);
dS = L/S;
rc = dS;

tic

% 1) Extend {x, f} to simplify periodic wrapping of boxes
[xe qe] = extend_one_layer(x, q, L, S);

% 2) Construct boxes over extended {x, f} 
[bx xe_bx_idx] = make_boxes(xe, -dS, L+dS, S+2);

% 3) Make neighbor list (note: indices into xe)
xe_nbh_list = make_nbh_list(xe_bx_idx, bx, N);

wtime.bx = toc;
tic

% 4) Prune neighbor list (Verlet)
if verlet
    xe_nbh_verlet = prune_list(xe, xe_nbh_list, rc);
else
    xe_nbh_verlet = xe_nbh_list;
end
wtime.verlet = toc;   

% 5) Compute interactions
tic
[nbh_vec nbh_len] = serialize_nbh_list(xe_nbh_verlet,'C');
phi = real_nbh_sum_mex(xe, qe, N, nbh_vec, nbh_len, xi);

wtime.sum = toc;

if checking
    % make sure everything conforms to assumptions
    check_bx_interior_26(bx);
    check_x_leading_in_y(x, xe);
    check_leading_x_interior_bx(xe_bx_idx, S+2, size(x, 1));

    if plotting
        plot3(xe(:,1), xe(:,2), xe(:,3), 'r*'), hold on, 
        plot3(x(:,1), x(:,2), x(:,3), '.'), 
        axis image
        bx_plot(1, bx, xe, -dS, L+dS, true);
        bx_plot(3, bx, xe, -dS, L+dS, true);
        
        nbh_plot(xe,xe_nbh_list,1,-dS, L+dS)
        nbh_plot(xe,xe_nbh_list,100,-dS, L+dS)

        nbh_plot(xe,xe_nbh_verlet,1,-dS, L+dS)
    end
end

if verbose
    print_info(x,xe,S,bx,xe_nbh_list,xe_nbh_verlet);
end

% ------------------------------------------------------------------------------
function pruned_list = prune_list(x,list,rc)
% remove from neighbor list the interactions farther than rc
N = length(list);
rc2 = rc^2;
for i=1:N
    p = x(i,:);
    q = x(list{i},:);
    d = (q(:,1)-p(1)).^2 + (q(:,2)-p(2)).^2 + (q(:,3)-p(3)).^2;
    pruned_list{i} = list{i}( d<= rc2 );
end

% ------------------------------------------------------------------------------
function x_nbh_list = make_nbh_list(x_bx_idx, bx, N)

for j=1:N
    
    % own box
    own_bx = x_bx_idx(j);
    
    % concatenate all negibors
    q = [bx(own_bx).x_idx; cat(1,bx( bx(own_bx).nbh(1:26) ).x_idx)];

    % remove self
    %q = setxor(q,j);
    
    % sort
    %q = sort(q);

    % done
    x_nbh_list{j} = q;
end

% ------------------------------------------------------------------------------
function [bx x_bx_idx] = make_boxes(x, a, b, S)
dS = (b-a)/S;
x_bx_sub = ceil((x-a)/dS);
x_bx_idx = sub2ind([S S S], x_bx_sub(:,1), x_bx_sub(:,2), x_bx_sub(:,3));

% construct box subscripts (needed???)
[sx sy sz] = ndgrid(1:S, 1:S, 1:S);
bx_sub = [sx(:) sy(:) sz(:)];

% neighbour mask
[nx ny nz] = ndgrid([-1 0 1], [-1 0 1], [-1 0 1]);
nbh = [nx(:) ny(:) nz(:)];
nbh( sum(abs(nbh), 2) < eps, :) = []; % self
assert(size(nbh,1) == 26) % have 26 neighbors 

% populate boxes
bx=repmat(struct('x_idx',[],'coords',[],'nbh',[]),[S^3,1]);
for i = 1:S^3

    % indices to particles that lie in my box
    bx(i).x_idx = int32(find(x_bx_idx == i));

    % current position in box decomposition (triple)
    bx(i).coords = bx_sub(i,:);

    % find neighbours (discard ones outside range [1, S])
    z = [nbh(:,1)+bx(i).coords(1) ...
         nbh(:,2)+bx(i).coords(2) ...
         nbh(:,3)+bx(i).coords(3)];
    keep = (z(:,1)>=1 & z(:,1) <=S) & ...
           (z(:,2)>=1 & z(:,2) <=S) & ...
           (z(:,3)>=1 & z(:,3) <=S);
    bx(i).nbh = sub2ind([S S S],z(keep,1),z(keep,2),z(keep,3));    
end

% ------------------------------------------------------------------------------
function [nbh_vec len] = serialize_nbh_list(nbh_list,s)
nbh_vec = int32(cat(1,nbh_list{:}));
switch s
    case 'C'
        nbh_vec = nbh_vec-1;
    case 'F'
        % do nothing
    otherwise
        error('must specify C or F(ortran) indexing')
end
len = int32(cellfun(@length,nbh_list));

% ------------------------------------------------------------------------------
function check_leading_x_interior_bx(x_bx_idx, S, N)
% check that all x lie in interior boxes
for j=1:N
    [ix iy iz] = ind2sub([S S S], x_bx_idx(j));
    assert(all( [ix iy iz] >= 2 ) && all( [ix iy iz] <= S-1 ))
end

function check_x_leading_in_y(x, y)
% check that x is the leading sub-array of y
for j = 1:size(x, 1)
    assert( all(abs( x(j,:)-y(j,:) ) < eps ) )
end

function check_bx_interior_26(bx)
% check that all interior boxes have 26 neighbors
Nbx = length(bx);
S = round(Nbx^(1/3));
for i = 2:S-1
    for j = 2:S-1
        for k = 2:S-1
            assert(length( bx(sub2ind([S S S], i, j, k) ).nbh ) == 26)
        end
    end
end

function bx_plot(idx, bx, x, a, b, plot_nbh)
% plot particles in a box and, possibly, all its neignbors
figure()
plot3(x(bx(idx).x_idx,1), x(bx(idx).x_idx,2), x(bx(idx).x_idx,3), '*'), 
axis([a b a b a b]), grid on, hold on
if plot_nbh
    for i = 1:length(bx(idx).nbh)
        plot3(x(bx( bx(idx).nbh(i) ).x_idx,1), ...
              x(bx( bx(idx).nbh(i) ).x_idx,2), ...
              x(bx( bx(idx).nbh(i) ).x_idx,3), 'r.'), 
    end
end

function nbh_plot(x,nbh,idx,a,b)
% plot x and all it's neighbors from list
figure()
plot3(x(idx,1),x(idx,2),x(idx,3),'r*')
axis([a b a b a b]), hold on
q=nbh{idx};
for j=1:length(q)
    plot3(x(q(j),1),x(q(j),2),x(q(j),3),'.')
end

function print_info(x,xe,S,bx,nbh_lst,nbh_verlet)
nb = numel(cat(1,nbh_lst{:}));
nv = numel(cat(1,nbh_verlet{:}));
N = size(x,1);
bxl = cellfun(@length,{bx.x_idx});
fprintf('[FAST RS] N %d \t Extended: %d\n',N,size(xe,1));
fprintf('          Boxes %d \t Extended: %d\n',S^3,length(bx));
fprintf('          Interactions(/N): BOX    %d \t%.1f\n', nb, nb/N);
fprintf('          Interactions(/N): PRUNE  %d \t%.1f\n', nv, nv/N);
fprintf('          Box occupancy: Max: %d   Min: %d   Sum: %d\n',...
    max(bxl),min(bxl),sum(bxl));

% ------------------------------------------------------------------------------
function [xe qe] = extend_one_layer(x, q, L, S)
% extend array x with one extra box layer to simplify periodicity

dS = L/S;
xe = x;
qe = q;

% x
ix1 = find(xe(:, 1) <= dS);
ix2 = find(xe(:, 1) >= L-dS);
x1 = [xe(ix1,1)+L xe(ix1,2) xe(ix1,3)];
x2 = [xe(ix2,1)-L xe(ix2,2) xe(ix2,3)];
xe = [xe; x1; x2];
qe = [qe; qe(ix1); qe(ix2)];

% y
ix1 = find(xe(:, 2) <= dS);
ix2 = find(xe(:, 2) >= L-dS);
x1 = [xe(ix1,1) xe(ix1,2)+L xe(ix1,3)];
x2 = [xe(ix2,1) xe(ix2,2)-L xe(ix2,3)];
xe = [xe; x1; x2];
qe = [qe; qe(ix1); qe(ix2)];

% z
ix1 = find(xe(:, 3) <= dS);
ix2 = find(xe(:, 3) >= L-dS);
x1 = [xe(ix1,1) xe(ix1,2) xe(ix1,3)+L];
x2 = [xe(ix2,1) xe(ix2,2) xe(ix2,3)-L];
xe = [xe; x1; x2];
qe = [qe; qe(ix1); qe(ix2)];
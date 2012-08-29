function [x f nvec] = generate_state(N,L, varargin)
%%Returning x: Nx3 matrix with coords inside box L, 
%%f: Nx3 matrix with point forces.
%%nvec: Nx3 matrix with associated normal vectors.
x = repmat(L,N,1).*rand(N,3);
f = 2*rand( N, 3) - 1;
nvec=rand( N, 3) - 0.5;
nabs=sqrt(nvec(:,1).^2+nvec(:,2).^2+nvec(:,3).^2); 
nvec=nvec./(nabs*ones(1,3));
%%Need to normalize normal vecs to length 1. 

if nargin==3 && strcmp(varargin{1},'neutral')
    % Charge neutral
    
    % Total sum from all but first 3 points
    S = zeros(3,3);
    for s=4:length(f);
        for l=1:3
            for m=1:3
                S(l,m) = S(l,m) + f(s,l)*nvec(s,m);
            end
        end
    end

    % Solve system for component l of first 4 points
    for l = 1:3
        A = [ nvec(1,:)' nvec(2,:)' nvec(3,:)' ];
        b = -[S(l,1); S(l,2); S(l,3)];
        y = A\b;
        f(1:3,l)=y;
    end
end

end
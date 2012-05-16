function [H1 H2 H3] = k_scaling(G,xi,box,eta)

M = size(G);
[H1 H2 H3] = deal(zeros(M(2:4)));
N = M(2); % Assuming cubic grid
if bitand(N,1)
    % odd
    kmax = (N-1)/2;
else
    kmax = N/2;
end
c = (1-eta)/(4*xi*xi);
for i2=0:N-1
    k(3) = -2.0*pi*(i2-kmax)/box(3);
    for i1=0:N-1
        k(2) = -2.0*pi*(i1-kmax)/box(2);
        for i0=0:N-1
            k(1) = -2.0*pi*(i0-kmax)/box(1);
            
            if(i0 ~= kmax || i1 ~= kmax || i2 ~= kmax)            
                k2 = k(1)*k(1)+k(2)*k(2)+k(3)*k(3);
                q = exp(-c*k2);

                B = 1i*stresslet_op_fd( k, xi)*q;
                block = G(1:9, i0+1, i1+1, i2+1);

                H1(i0+1, i1+1, i2+1) = B(1,:)*block;
                H2(i0+1, i1+1, i2+1) = B(2,:)*block;
                H3(i0+1, i1+1, i2+1) = B(3,:)*block;
            else
                % Continue
                disp('Skipping k=')
                disp(k);
            end        
        end
    end
end
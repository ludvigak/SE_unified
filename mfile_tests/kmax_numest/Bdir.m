function [ f ] = Bdir( j,l,m, theta, phi, k, r )

    khat = {cos(phi).*sin(theta)
            sin(phi).*sin(theta)
            cos(theta)};

    f = -pi*sin(k.*r.*cos(theta)) .* ...
        ((j==l)*khat{m} + (l==m)*khat{j} + (m==j)*khat{l} - 2*khat{j}.*khat{l}.*khat{m});
end


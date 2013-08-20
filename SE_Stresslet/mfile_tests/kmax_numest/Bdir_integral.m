function [ f ] = Bdir_integral( j,l,m,k,r)
    if j==3 && m==3 && l==3
        f = 4*pi^2*(k.*r.*(12+k.^2.*r.^2).*cos(k.*r)+3*(-4+k.^2.*r.^2).*sin(k.*r))./(k.^4.*r.^4);
    elseif ( (j==l)*(m==3) || (m==j)*(l==3) + (l==m)*(j==3) ) 
         f = 4*pi^2*(k.*r.*(-6+k.^2.*r.^2).*cos(k.*r)-3*(-2+k.^2.*r.^2).*sin(k.*r))./(k.^4.*r.^4);
    else
        f = 0;
    end
end


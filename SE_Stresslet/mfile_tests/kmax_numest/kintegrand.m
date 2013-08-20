function [ f ] = kintegrand(j,l,m,k,r,xi)
    % full k integrand, including volume k^2
    if r==0
        f = 0.*k.*r;
    else
        f =  Bdir_integral( j,l,m,k,r).*...
            exp(-k.^2/4/xi^2).*(8+2*k.^2/xi^2+k.^4/xi^4).*k;
    end
end


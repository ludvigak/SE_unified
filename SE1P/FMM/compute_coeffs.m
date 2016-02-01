function a = compute_coeffs(P, z,scale,center)
%% FIXME: This should be vectroized for particles
b = zeros(P+1,P+1);
a = zeros(P+1,P+1);
c = a;

gamma_ = 0.577215664901532860606512090082402431042;
x = (z(1)-center(1))*scale
y = (z(2)-center(2))*scale

x0 = x;
y0 = y;
x2y2 = x.^2+y.^2;
%  indices are 1 based

if(x2y2>100000)
    x2y2 = x.^2+y.^2;
    %  indices are 1 based

    b(1,1) = x2y2;
    b(2,1) = 2*x;
    b(1,2) = 2*y;
    b(2,2) = 0;
    b(3,1) = 1;
    b(1,3) = 1;
    %a(1,1) = expint(x2y2) + log(x^2+y^2) + gamma_;
    a(1,1) = log(x^2+y^2);
    a(2,1) = 2*x/x2y2;
    a(1,2) = 2*y/x2y2;
    a(2,2) = -4*x*y/(x2y2)^2;
    
    for k1 = 2:P
        k2 = 0; s = k1+k2;
        a(k1+1,k2+1) = (-2*(s-1)*x*a(k1,1) - (s-2)*a(k1-1,1) + b(k1+1,1)*s)/s/x2y2;
        a(k2+1,k1+1) = (-2*(s-1)*y*a(1,k1) - (s-2)*a(1,k1-1) + b(1,k1+1)*s)/s/x2y2;
        
        k2 = 1; s = k1+k2;     
        
        a(k1+1,k2+1) = (b(k1+1,k2+1)-2*x*a(k1,k2+1)-a(k1-1,k2+1))/x2y2;
        a(k2+1,k1+1) = (b(k1+1,k2+1)-2*y*a(k2+1,k1)-a(k2+1,k1-1))/x2y2;
        a(k1+1,k2+1) = (-2*(s-1)* (x*a(k1,k2+1) + y*a(k1+1,k2)) - (s-2)*a(k1-1,k2+1) + b(k1+1,k2+1)*s)/s/x2y2;
        a(k2+1,k1+1) = (-2*(s-1)* (x*a(k2,k1+1) + y*a(k2+1,k1)) - (s-2)*a(k2+1,k1-1) + b(k2+1,k1+1)*s)/s/x2y2;

        for k2 =2:P
            if( (k1+k2)>P)
                continue;
            end
            s = k1+k2;
            a(k1+1,k2+1) = (-2*(s-1) * (x*a(k1  ,k2+1) + y*a(k1+1,k2  )) ...
                -  (s-2) * (  a(k1-1,k2+1) +   a(k1+1,k2-1)) ...
                + s    *    b(k1+1,k2+1)                 )/s/x2y2;
            a(k2+1,k1+1) = (-2*(s-1) * (x*a(k2  ,k1+1) + y*a(k2+1,k1  )) ...
                -  (s-2) * (  a(k2-1,k1+1) +   a(k2+1,k1-1)) ...
                + s    *    b(k2+1,k1+1)                 )/s/x2y2;
            
        end
    end

else
    %% construct bk with recursion
    %  indices are 1 based
    exp2 = exp(-x2y2);

    % divison by 0! and 1! are ignored
    b(1,1) = exp2;
    b(2,1) = -2*x*exp2;
    b(1,2) = -2*y*exp2;
    b(2,2) = 4*x*y*exp2;

    for k1 = 2:P
        b(k1+1,1) = -2/k1 *( x*b(k1,1) + b(k1-1,1) );
        b(k1+1,2) = -2/k1 *( x*b(k1,2) + b(k1-1,2) );
        b(1,k1+1) = -2/k1 *( y*b(1,k1) + b(1,k1-1) );
        b(2,k1+1) = -2/k1 *( y*b(2,k1) + b(2,k1-1) );
        for k2 = 2:P
            b(k1+1,k2+1) = -2/k2 * ( y*b(k1+1,k2)+b(k1+1,k2-1) );
        end
    end

    b(1,1) = b(1,1) + x2y2;
    b(2,1) = b(2,1) + 2*x;
    b(1,2) = b(1,2) + 2*y;
    b(2,2) = b(2,2) +  0;
    b(3,1) = b(3,1) + 1;
    b(1,3) = b(1,3) + 1;

    %% construct ak with recursion
    % indices are 1 based
    % Adding gamma to this function wont affect the performance.
    a(1,1) = expint(x2y2)+log(x^2+y^2) + gamma_;
    a(2,1) = 2*x*(1-exp2)/x2y2;
    a(1,2) = 2*y*(1-exp2)/x2y2;
    a(2,2) = 4*x*y*(exp2*(x2y2+1) - 1)/(x2y2)^2;
    
    for k1 = 2:P
        k2 = 0; s = k1+k2;
        a(k1+1,k2+1) = (-2*(s-1)*x*a(k1,1) - (s-2)*a(k1-1,1) + b(k1+1,1)*s)/s/x2y2;
        a(k2+1,k1+1) = (-2*(s-1)*y*a(1,k1) - (s-2)*a(1,k1-1) + b(1,k1+1)*s)/s/x2y2;

    %    a(k1+1,k2+1) = (b(k1+1,1)-2*x*(1-1/k1)*a(k1,1)-(1-2/k1)*a(k1-1,1))/x2y2;
    %    a(k2+1,k1+1) = (b(1,k1+1)-2*y*(1-1/k1)*a(1,k1)-(1-2/k1)*a(1,k1-1))/x2y2;

        k2 = 1; s = k1+k2;
         a(k1+1,k2+1) = (-2*(s-1)* (x*a(k1,k2+1) + y*a(k1+1,k2)) - (s-2)*a(k1-1,k2+1) + b(k1+1,k2+1)*s)/s/x2y2;
         a(k2+1,k1+1) = (-2*(s-1)* (x*a(k2,k1+1) + y*a(k2+1,k1)) - (s-2)*a(k2+1,k1-1) + b(k2+1,k1+1)*s)/s/x2y2;

         %a(k1+1,k2+1) = (b(k1+1,k2+1)-2*x*a(k1,k2+1)-a(k1-1,k2+1))/x2y2;
         %a(k2+1,k1+1) = (b(k1+1,k2+1)-2*y*a(k2+1,k1)-a(k2+1,k1-1))/x2y2;
        for k2 =2:P
            s = k1+k2;
            a(k1+1,k2+1) = (-2*(s-1) * (x*a(k1  ,k2+1) + y*a(k1+1,k2  )) ...
                -  (s-2) * (  a(k1-1,k2+1) +   a(k1+1,k2-1)) ...
                + s    *    b(k1+1,k2+1)                 )/s/x2y2;
            a(k2+1,k1+1) = (-2*(s-1) * (x*a(k2  ,k1+1) + y*a(k2+1,k1  )) ...
                -  (s-2) * (  a(k2-1,k1+1) +   a(k2+1,k1-1)) ...
                + s    *    b(k2+1,k1+1)                 )/s/x2y2;
        end
    end
    
    
end

test_coeff = 0;
if (test_coeff == 1)
    syms x y
    f = expint(x^2+y^2)+log(x^2+y^2)+gamma_;
    f = log(x^2+y^2);
    g = exp(-x^2-y^2);
    
    k1 = 2;
    k2 = 0;
    kfact = factorial(k1)*factorial(k2);
    fder = diff(diff(f,'x',k1),'y',k2);
    gder = diff(diff(g,'x',k1),'y',k2);
    x = x0;
    y = y0;
    fprintf('b:exact\t\tcoeff\t\tdiff\n');
    der = eval(gder);
    fprintf('  %f\t%f\t%g\n',der/kfact,b(k1+1,k2+1),der/kfact - b(k1+1,k2+1));
    fprintf('a:exact\t\tcoeff\t\tdiff\n');
    der = eval(fder);
    fprintf('  %.16f\t%.16f\t%g\n',der/kfact,a(k1+1,k2+1),der/kfact-a(k1+1,k2+1));
end
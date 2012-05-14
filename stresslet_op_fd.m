function Bi = stresslet_op_fd( k, xi)
% Returns imaginary part of Bi(k, xi)
% ie B = 1i*Bi;
k2 = k(1)^2+k(2)^2+k(3)^2;
k2i = 1/k2;
c = k2/xi^2;

psi=-pi*(8.0+2*c+c^2)*k2i;

Bi = zeros(3,3,3);

% Explicit code
% for m=1:3
%     for l=1:3
%         for j=1:3
%             Bi(j,l,m) = psi*(...
%                 (j==l)*k(m)+...
%                 (l==m)*k(j)+...
%                 (m==j)*k(l)-...
%                 2*k(j)*k(l)*k(m)*k2i);
%         end
%     end
% end

% Unrolled
Bi(1,1,1)=psi*(k(1)+k(1)+k(1)-2*k(1)*k(1)*k(1)*k2i);
Bi(2,1,1)=psi*(k(2)-2*k(2)*k(1)*k(1)*k2i);
Bi(3,1,1)=psi*(k(3)-2*k(3)*k(1)*k(1)*k2i);
Bi(1,2,1)=psi*(k(2)-2*k(1)*k(2)*k(1)*k2i);
Bi(2,2,1)=psi*(k(1)-2*k(2)*k(2)*k(1)*k2i);
Bi(3,2,1)=psi*(-2*k(3)*k(2)*k(1)*k2i);
Bi(1,3,1)=psi*(k(3)-2*k(1)*k(3)*k(1)*k2i);
Bi(2,3,1)=psi*(-2*k(2)*k(3)*k(1)*k2i);
Bi(3,3,1)=psi*(k(1)-2*k(3)*k(3)*k(1)*k2i);
Bi(1,1,2)=psi*(k(2)-2*k(1)*k(1)*k(2)*k2i);
Bi(2,1,2)=psi*(k(1)-2*k(2)*k(1)*k(2)*k2i);
Bi(3,1,2)=psi*(-2*k(3)*k(1)*k(2)*k2i);
Bi(1,2,2)=psi*(k(1)-2*k(1)*k(2)*k(2)*k2i);
Bi(2,2,2)=psi*(k(2)+k(2)+k(2)-2*k(2)*k(2)*k(2)*k2i);
Bi(3,2,2)=psi*(k(3)-2*k(3)*k(2)*k(2)*k2i);
Bi(1,3,2)=psi*(-2*k(1)*k(3)*k(2)*k2i);
Bi(2,3,2)=psi*(k(3)-2*k(2)*k(3)*k(2)*k2i);
Bi(3,3,2)=psi*(k(2)-2*k(3)*k(3)*k(2)*k2i);
Bi(1,1,3)=psi*(k(3)-2*k(1)*k(1)*k(3)*k2i);
Bi(2,1,3)=psi*(-2*k(2)*k(1)*k(3)*k2i);
Bi(3,1,3)=psi*(k(1)-2*k(3)*k(1)*k(3)*k2i);
Bi(1,2,3)=psi*(-2*k(1)*k(2)*k(3)*k2i);
Bi(2,2,3)=psi*(k(3)-2*k(2)*k(2)*k(3)*k2i);
Bi(3,2,3)=psi*(k(2)-2*k(3)*k(2)*k(3)*k2i);
Bi(1,3,3)=psi*(k(1)-2*k(1)*k(3)*k(3)*k2i);
Bi(2,3,3)=psi*(k(2)-2*k(2)*k(3)*k(3)*k2i);
Bi(3,3,3)=psi*(k(3)+k(3)+k(3)-2*k(3)*k(3)*k(3)*k2i);

% Unrolled and simplified
% Bi(1,1,1)=psi*(k(1)+k(1)+k(1)-2*k(1)*k(1)*k(1)*k2i);
% Bi(2,1,1)=psi*(k(2)-2*k(2)*k(1)*k(1)*k2i);
% Bi(3,1,1)=psi*(k(3)-2*k(3)*k(1)*k(1)*k2i);
% Bi(1,2,1)=Bi(2,1,1);
% Bi(2,2,1)=psi*(k(1)-2*k(2)*k(2)*k(1)*k2i);
% Bi(3,2,1)=psi*(-2*k(3)*k(2)*k(1)*k2i);
% Bi(1,3,1)=Bi(3,1,1);
% Bi(2,3,1)=Bi(3,2,1);
% Bi(3,3,1)=psi*(k(1)-2*k(3)*k(3)*k(1)*k2i);
% Bi(1,1,2)=Bi(1,2,1);
% Bi(2,1,2)=Bi(2,2,1);
% Bi(3,1,2)=Bi(3,2,1);
% Bi(1,2,2)=Bi(2,1,2);
% Bi(2,2,2)=psi*(k(2)+k(2)+k(2)-2*k(2)*k(2)*k(2)*k2i);
% Bi(3,2,2)=psi*(k(3)-2*k(3)*k(2)*k(2)*k2i);
% Bi(1,3,2)=Bi(3,1,2);
% Bi(2,3,2)=Bi(3,2,2);
% Bi(3,3,2)=psi*(k(2)-2*k(3)*k(3)*k(2)*k2i);
% Bi(1,1,3)=Bi(1,3,1);
% Bi(2,1,3)=Bi(2,3,1);
% Bi(3,1,3)=Bi(3,3,1);
% Bi(1,2,3)=Bi(2,1,3);
% Bi(2,2,3)=Bi(2,3,2);
% Bi(3,2,3)=Bi(3,3,2);
% Bi(1,3,3)=Bi(3,1,3);
% Bi(2,3,3)=Bi(3,2,3);
% Bi(3,3,3)=psi*(k(3)+k(3)+k(3)-2*k(3)*k(3)*k(3)*k2i);


% % Shorter way, slower (!)
% Bi(:,:,1) = eye(3)-2*(khat')*khat;
% khat = khat*psi;
% for j=3:-1:1;
%     Bi(:,:,j) = khat(j)*Bi(:,:,1);
%     Bi(j,:,j) = Bi(j,:,j) + khat;
%     Bi(:,j,j) = Bi(:,j,j) + khat';
% end

end


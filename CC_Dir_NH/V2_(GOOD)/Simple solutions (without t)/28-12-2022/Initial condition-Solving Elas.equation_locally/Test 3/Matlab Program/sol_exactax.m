function sol=sol_exactax(xx,yy,tt,comp)

global perm

%       u(xx,yy,tt) = {u1(xx,yy,tt),u2(xx,yy,tt)}

%         p(xx,yy,tt)

%         z(xx,yy,tt)= {z1(xx,yy,tt),z2xx,yy,tt)}

%         gamma = [0 r;-r 0]
%               = r

 if comp==1         % Soluci�n anal�tica de u1(xx,yy,tt)
     sol=(1 - xx)*xx*yy;    
 elseif comp==2     % Soluci�n anal�tica de u2(xx,yy,tt)
     sol=xx*(1 - yy)*yy;
 elseif comp==3     % Soluci�n anal�tica de p(xx,yy,tt)
     sol=xx*yy*(1 - xx*yy);
 elseif comp==4     % Soluci�n anal�tica de z1(xx,yy,tt)
     sol=perm*yy*(-1 + 2*xx*yy);
 elseif comp==5     % Soluci�n anal�tica de z2(xx,yy,tt)
     sol=perm*xx*(-1 + 2*xx*yy);
 else               % Soluci�n anal�tica de gamma (1st row & 2nd column)
     sol=(xx - xx^2 + (-1 + yy)*yy)/2.;
 end

return
end
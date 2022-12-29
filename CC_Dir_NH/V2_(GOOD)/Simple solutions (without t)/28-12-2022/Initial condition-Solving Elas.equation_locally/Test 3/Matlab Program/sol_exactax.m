function sol=sol_exactax(xx,yy,tt,comp)

global perm

%       u(xx,yy,tt) = {u1(xx,yy,tt),u2(xx,yy,tt)}

%         p(xx,yy,tt)

%         z(xx,yy,tt)= {z1(xx,yy,tt),z2xx,yy,tt)}

%         gamma = [0 r;-r 0]
%               = r

 if comp==1         % Solución analítica de u1(xx,yy,tt)
     sol=(1 - xx)*xx*yy;    
 elseif comp==2     % Solución analítica de u2(xx,yy,tt)
     sol=xx*(1 - yy)*yy;
 elseif comp==3     % Solución analítica de p(xx,yy,tt)
     sol=xx*yy*(1 - xx*yy);
 elseif comp==4     % Solución analítica de z1(xx,yy,tt)
     sol=perm*yy*(-1 + 2*xx*yy);
 elseif comp==5     % Solución analítica de z2(xx,yy,tt)
     sol=perm*xx*(-1 + 2*xx*yy);
 else               % Solución analítica de gamma (1st row & 2nd column)
     sol=(xx - xx^2 + (-1 + yy)*yy)/2.;
 end

return
end
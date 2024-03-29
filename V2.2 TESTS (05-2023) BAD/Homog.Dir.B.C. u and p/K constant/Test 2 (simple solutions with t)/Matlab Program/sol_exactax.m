function sol=sol_exactax(xx,yy,tt,comp)

global perm

%       u(xx,yy,tt) = {u1(xx,yy,tt),u2(xx,yy,tt)}

%         p(xx,yy,tt)

%         z(xx,yy,tt)= {z1(xx,yy,tt),z2xx,yy,tt)}

%         gamma = [0 r;-r 0]
%               = r

 if comp==1         % Solución analítica de u1(xx,yy,tt)
     sol=tt*(1 - xx)*xx*(1 - yy)*yy;
 elseif comp==2     % Solución analítica de u2(xx,yy,tt)
     sol=tt*(1 - xx)*xx*(1 - yy)*yy;
 elseif comp==3     % Solución analítica de p(xx,yy,tt)
%      sol=tt*(1 - xx)*xx*(1 - yy)*yy;
     sol=tt.*(1 - xx).*xx.*(1 - yy).*yy;
 elseif comp==4     % Solución analítica de z1(xx,yy,tt)
     sol=-(perm*tt*(-1 + 2*xx)*(-1 + yy)*yy);
 elseif comp==5     % Solución analítica de z2(xx,yy,tt)
     sol=-(perm*tt*(-1 + xx)*xx*(-1 + 2*yy));
 else               % Solución analítica de gamma (1st row & 2nd column)
     sol=(tt*(xx + (-1 + yy)*yy - 2*xx*yy^2 + xx^2*(-1 + 2*yy)))/2.;
 end

return
end
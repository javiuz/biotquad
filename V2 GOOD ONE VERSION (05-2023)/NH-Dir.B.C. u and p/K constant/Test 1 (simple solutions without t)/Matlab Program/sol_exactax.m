function sol=sol_exactax(xx,yy,tt,comp)

global perm

%         p(xx,yy,tt)
%       u(xx,yy,tt) = {u1(xx,yy,tt),u2(xx,yy,tt)}


%         z(xx,yy,tt)= {z1(xx,yy,tt),z2xx,yy,tt)}

%         gamma = [0 r;-r 0]
%               = r

 if comp==1         % Solución analítica de u1(xx,yy,tt)
     sol=xx;    
 elseif comp==2     % Solución analítica de u2(xx,yy,tt)
     sol=yy;
 elseif comp==3     % Solución analítica de p(xx,yy,tt)
     sol=yy;
 elseif comp==4     % Solución analítica de z1(xx,yy,tt)
     sol=0;
 elseif comp==5     % Solución analítica de z2(xx,yy,tt)
     sol=-perm;
 else               % Solución analítica de gamma (1st row & 2nd column)
     sol=0;
 end

return
end
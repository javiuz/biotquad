function sol=sol_exactax(xx,yy,tt,comp)

global perm

%       u(xx,yy,tt) = {u1(xx,yy,tt),u2(xx,yy,tt)}

%         p(xx,yy,tt)

%         z(xx,yy,tt)= {z1(xx,yy,tt),z2xx,yy,tt)}

%         gamma = [0 r;-r 0]
%               = r

 if comp==1         % Solución analítica de u1(xx,yy,tt)
     sol=(xx^2 + xx^3*yy^4 + cos(1 - yy)*sin((1 - xx)*(1 - yy)));
 elseif comp==2     % Solución analítica de u2(xx,yy,tt)
     sol=((1 - yy)^2 + (1 - xx)^4*(1 - yy)^3 + cos(xx*yy)*sin(xx));
 elseif comp==3     % Solución analítica de p(xx,yy,tt)
     sol=(10 + cos(pi*yy)*sin(pi*xx));
 elseif comp==4     % Solución analítica de z1(xx,yy,tt)
     sol=-(perm*pi*cos(pi*xx)*cos(pi*yy));
 elseif comp==5     % Solución analítica de z2(xx,yy,tt)
     sol=perm*pi*sin(pi*xx)*sin(pi*yy);
 else               % Solución analítica de gamma (1st row & 2nd column)
     sol=((4*(-1 + xx)^3*(-1 + yy)^3 + 4*xx^3*yy^3 + (-1 + xx)*cos(1 - yy)*cos((-1 + xx)*(-1 + yy)) - ...
         cos(xx)*cos(xx*yy) + sin(1 - yy)*sin((-1 + xx)*(-1 + yy)) + yy*sin(xx)*sin(xx*yy)))/2.;
 end

return
end
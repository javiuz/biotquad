function sol=sol_exactax(xx,yy,tt,comp)

global perm

%       u(xx,yy) = {x,y}

%         u1(xx) = xx
%         u2(yy) = yy

%         p(xx,tt) = tt*xx

%         z(tt)= {z1(tt),z2}
%           z1(tt)=-perm*tt
%           z2=0

%         gamma = [0 r;-r 0]
%           r=0

 if comp==1         % Solución analítica de u1(xx)
     sol=xx*(1-xx)*yy*(1-yy);    
 elseif comp==2     % Solución analítica de u2(yy)
     sol=xx*(xx-1)*yy*(1-yy);
 elseif comp==3     % Solución analítica de p(xx,tt)
     sol=xx*(1-xx)*yy*(1-yy);
 elseif comp==4     % Solución analítica de z1(tt)
     sol=-perm*(-1 + 2*xx)*(-1 + yy)*yy;
 elseif comp==5     % Solución analítica de z2
     sol=-perm*(-1 + xx)*xx*(-1 + 2 *yy);
 else               % Solución analítica de gamma (1st row & 2nd column)
     sol=(yy + 2*xx*(-1 + yy)*yy - yy^2 + (-1 + xx)*xx*(-1 + 2*yy))/2.;
 end

return
end
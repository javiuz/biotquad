function sol=sol_exactax(xx,yy,tt,comp)

global perm

%         p(xx,yy,tt)
%       u(xx,yy,tt) = {u1(xx,yy,tt),u2(xx,yy,tt)}


%         z(xx,yy,tt)= {z1(xx,yy,tt),z2xx,yy,tt)}

%         gamma = [0 r;-r 0]
%               = r

 if comp==1         % Soluci�n anal�tica de u1(xx,yy,tt)
     sol=xx;    
 elseif comp==2     % Soluci�n anal�tica de u2(xx,yy,tt)
     sol=yy;
 elseif comp==3     % Soluci�n anal�tica de p(xx,yy,tt)
     sol=yy;
 elseif comp==4     % Soluci�n anal�tica de z1(xx,yy,tt)
     sol=0;
 elseif comp==5     % Soluci�n anal�tica de z2(xx,yy,tt)
     sol=-perm;
 else               % Soluci�n anal�tica de gamma (1st row & 2nd column)
     sol=0;
 end

return
end
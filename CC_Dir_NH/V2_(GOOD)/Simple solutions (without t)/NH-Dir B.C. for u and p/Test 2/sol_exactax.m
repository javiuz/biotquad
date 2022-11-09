function sol=sol_exactax(xx,yy,~,comp)

global perm

%       u(xx,yy) = {x,y}

%         u1(xx) = xx
%         u2(yy) = yy

%         p(xx,yy) = yy*xx

%         z(tt)= [z1,z2]'
%           z1(tt)=-perm*yy
%           z2=-perm*xx

%         gamma = [0 r;-r 0]
%           r=0

 if comp==1         % Soluci�n anal�tica de u1(xx)
     sol=xx;    
 elseif comp==2     % Soluci�n anal�tica de u2(yy)
     sol=yy;
 elseif comp==3     % Soluci�n anal�tica de p(xx,tt)
     sol=yy*xx;
 elseif comp==4     % Soluci�n anal�tica de z1(tt)
     sol=-perm*yy;
 elseif comp==5     % Soluci�n anal�tica de z2
     sol=-perm*xx;
 else               % Soluci�n anal�tica de gamma (1st row & 2nd column)
     sol=0;
 end

return
end
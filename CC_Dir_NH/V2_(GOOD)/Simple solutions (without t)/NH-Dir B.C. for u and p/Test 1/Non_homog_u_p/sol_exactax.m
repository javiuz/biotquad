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
     sol=xx; 
%      sol=0; 
 elseif comp==2     % Solución analítica de u2(yy)
     sol=yy;
%      sol=0;
 elseif comp==3     % Solución analítica de p(xx,tt)
     sol=xx;
%      sol=yy;
 elseif comp==4     % Solución analítica de z1(tt)
     sol=-perm;
%      sol=0;
 elseif comp==5     % Solución analítica de z2
     sol=0;
 else               % Solución analítica de gamma (1st row & 2nd column)
     sol=0;
 end

return
end
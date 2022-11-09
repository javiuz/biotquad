function sol=sol_exactax_sigma(xx,yy,tt,~,row,col)

global alpha lambda mu

%       sigma(xx,yy,tt) = [sigma(1,1) sigma(1,2);sigma(2,1) sigma(2,2)]

 if row==1 && col==1       % Solución analítica de sigma (1,1) (xx,yy,tt)
     sol=alpha*xx*yy*(-1 + xx + yy - xx*yy) + 2*(lambda + mu);     
 elseif row==1 && col==2   % Solución analítica de sigma (1,2) (xx,yy,tt)
     sol=0;
 elseif row==2 && col==1   % Solución analítica de sigma (2,1) (xx,yy,tt)
     sol=0;
 else                       % Solución analítica de sigma (2,2) (xx,yy,tt)
     sol=alpha*xx*yy*(-1 + xx + yy - xx*yy) + 2*(lambda + mu);    
 end

return
end
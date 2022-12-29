function sol=sol_exactax_sigma(xx,yy,tt,~,row,col)

global alpha lambda mu

%       sigma(xx,yy,tt) = [sigma(1,1) sigma(1,2);sigma(2,1) sigma(2,2)]

 if row==1 && col==1       % Solución analítica de sigma (1,1) (xx,yy,tt)
     sol=mu*(2 - 4*xx)*yy + lambda*(xx + yy - 4*xx*yy) + alpha*xx*yy*(-1 + xx*yy);     
 elseif row==1 && col==2   % Solución analítica de sigma (1,2) (xx,yy,tt)
     sol=mu*(xx - xx^2 + yy - yy^2);
 elseif row==2 && col==1   % Solución analítica de sigma (2,1) (xx,yy,tt)
     sol=mu*(xx - xx^2 + yy - yy^2);
 else                       % Solución analítica de sigma (2,2) (xx,yy,tt)
     sol=mu*xx*(2 - 4*yy) + lambda*(xx + yy - 4*xx*yy) + alpha*xx*yy*(-1 + xx*yy);   
 end

return
end
function sol=sol_exactax_sigma(xx,yy,tt,~,row,col)

global alpha lambda mu

%       sigma(xx,yy,tt) = [sigma(1,1) sigma(1,2);sigma(2,1) sigma(2,2)]

 if row==1 && col==1       % Solución analítica de sigma (1,1) (xx,yy,tt)
     sol=(-(alpha*(-1 + xx)*xx) + mu*(-2 + 4*xx))*(-1 + yy)*yy + lambda*(-xx + xx^2 + yy - 2*xx^2*yy - yy^2 + 2*xx*yy^2);     
 elseif row==1 && col==2   % Solución analítica de sigma (1,2) (xx,yy,tt)
     sol=mu*(xx + (-1 + yy)*yy - 2*xx*yy^2 + xx^2*(-1 + 2*yy));
 elseif row==2 && col==1   % Solución analítica de sigma (2,1) (xx,yy,tt)
     sol=mu*(xx + (-1 + yy)*yy - 2*xx*yy^2 + xx^2*(-1 + 2*yy));
 else                       % Solución analítica de sigma (2,2) (xx,yy,tt)
     sol=lambda*(-xx + xx^2 + yy - 2*xx^2*yy - yy^2 + 2*xx*yy^2) - (-1 + xx)*xx*(alpha*(-1 + yy)*yy + mu*(-2 + 4*yy));    
 end

return
end
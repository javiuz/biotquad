function sol=sol_exactax_sigma(xx,yy,tt,~,row,col)

global alpha lambda mu

%       sigma(xx,yy,tt) = [sigma(1,1) sigma(1,2);sigma(2,1) sigma(2,2)]

 if row==1 && col==1       % Solución analítica de sigma (1,1) (xx,yy,tt)
     sol=lambda*(-1 + 2*xx)*(-1 + yy)*yy + 2*mu*(-1 + 2*xx)*(-1 + yy)*yy + lambda*(-1 + xx)*xx*(-1 + 2*yy);     
 elseif row==1 && col==2   % Solución analítica de sigma (1,2) (xx,yy,tt)
     sol=mu*(-((-1 + yy)*yy) + xx^2*(-1 + 2*yy) + xx*(1 - 4*yy + 2*yy^2));
 elseif row==2 && col==1   % Solución analítica de sigma (2,1) (xx,yy,tt)
     sol=mu*(-((-1 + yy)*yy) + xx^2*(-1 + 2*yy) + xx*(1 - 4*yy + 2*yy^2));
 else                       % Solución analítica de sigma (2,2) (xx,yy,tt)
     sol=lambda*(-1 + 2*xx)*(-1 + yy)*yy + lambda*(-1 + xx)*xx*(-1 + 2*yy) + 2*mu*(-1 + xx)*xx*(-1 + 2*yy);   
 end

return
end
function sol=sol_exactax_sigma(xx,yy,tt,nu,row,col)

global alpha lambda mu

%       sigma(xx,yy,tt) = [sigma(1,1) sigma(1,2);sigma(2,1) sigma(2,2)]

 if row==1 && col==1       % Solución analítica de sigma (1,1) (xx,yy,tt)
     sol=(2*mu*(xx*(2 + 3*xx*yy^4) + (-1 + yy)*cos(1 - yy)*cos((-1 + xx)*(-1 + yy))) - ...
         alpha*(10 + cos(pi*yy)*sin(pi*xx)) + lambda*...
         (2*(-1 + yy) - 3*(-1 + xx)^4*(-1 + yy)^2 + xx*(2 + 3*xx*yy^4) + ...
         (-1 + yy)*cos(1 - yy)*cos((-1 + xx)*(-1 + yy)) - xx*sin(xx)*sin(xx*yy)));     
 elseif row==1 && col==2   % Solución analítica de sigma (1,2) (xx,yy,tt)
     sol=mu*(-4*(-1 + xx)^3*(-1 + yy)^3 + 4*xx^3*yy^3 + (-1 + xx)*cos(1 - yy)*cos((-1 + xx)*(-1 + yy)) + ...
         cos(xx)*cos(xx*yy) + sin(1 - yy)*sin((-1 + xx)*(-1 + yy)) - yy*sin(xx)*sin(xx*yy));
 elseif row==2 && col==1   % Solución analítica de sigma (2,1) (xx,yy,tt)
     sol=mu*(-4*(-1 + xx)^3*(-1 + yy)^3 + 4*xx^3*yy^3 + (-1 + xx)*cos(1 - yy)*cos((-1 + xx)*(-1 + yy)) + ...
         cos(xx)*cos(xx*yy) + sin(1 - yy)*sin((-1 + xx)*(-1 + yy)) - yy*sin(xx)*sin(xx*yy));
 else                       % Solución analítica de sigma (2,2) (xx,yy,tt)
     sol=(-(alpha*(10 + cos(pi*yy)*sin(pi*xx))) + ...
         2*mu*(2*(-1 + yy) - 3*(-1 + xx)^4*(-1 + yy)^2 - xx*sin(xx)*sin(xx*yy)) + ...
         lambda*(2*(-1 + yy) - 3*(-1 + xx)^4*(-1 + yy)^2 + xx*(2 + 3*xx*yy^4) + ...
         (-1 + yy)*cos(1 - yy)*cos((-1 + xx)*(-1 + yy)) - xx*sin(xx)*sin(xx*yy)));   
 end

return
end
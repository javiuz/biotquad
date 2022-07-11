function sol=sol_exactax(xx,yy,tt,comp)

%       u(xx,yy,tt) = {e^t (x^2 + x^3 y^4 + Cos[1 - y] Sin[(1 - x) (1 - y)]), 
%                     e^t ((1 - y)^2 + (1 - x)^4 (1 - y)^3 + Cos[x y] Sin[x])}

%         u1(xx,yy,tt) = e^t (x^2 + x^3 y^4 + Cos[1 - y] Sin[(1 - x) (1 - y)])
%         u2(xx,yy,tt) = e^t ((1 - y)^2 + (1 - x)^4 (1 - y)^3 + Cos[x y] Sin[x])

%         p(xx,yy,tt) = e^t (10 + Cos[pi y] Sin[pi x])

%         z(xx,yy,tt)= {z1(xx,yy,tt),z2(xx,yy,tt)}
%           z1(xx,yy,tt)=-(E**tt*pi*((1 + xx)**2 + yy**2)*cos(pi*xx)*cos(pi*yy)) + E**tt*pi*sin(pi*xx)*sin(pi*yy)*sin(xx*yy)
%           z2(xx,yy,tt)=E**tt*pi*(1 + xx)**2*sin(pi*xx)*sin(pi*yy) - E**tt*pi*cos(pi*xx)*cos(pi*yy)*sin(xx*yy)

%         gamma(xx,yy,tt) = [0 r(xx,yy,tt);-r(xx,yy,tt) 0]
%           r(xx,yy,tt)=1/2 (e^t (4 x^3 y^3 + (-1 + x) Cos[1 - y] Cos[(1 - x) (1 - y)] + Sin[1 - y] Sin[(1 - x) (1 - y)]) - e^t (-4 (1 - x)^3 (1 - y)^3 + Cos[x] Cos[x y] - y Sin[x] Sin[x y]))

 if comp==1         % Solución analítica de u1(xx,yy,tt)
     sol=exp(tt)*(xx^2 + xx^3*yy^4 + cos(1-yy)*sin((1-xx)*(1-yy)));
     
 elseif comp==2     % Solución analítica de u2(xx,yy,tt)
     sol=exp(tt)*((1-yy)^2 + (1-xx)^4*(1-yy)^3 + cos(xx*yy)*sin(xx));
 elseif comp==3     % Solución analítica de p(xx,yy,tt)
     sol=exp(tt)*(10 + cos(pi*yy)*sin(pi*xx));
 elseif comp==4     % Solución analítica de z1(xx,yy,tt)
     sol=-(exp(tt)*pi*((1 + xx)^2 + yy^2)*cos(pi*xx)*cos(pi*yy)) + exp(tt)*pi*sin(pi*xx)*sin(pi*yy)*sin(xx*yy);
 elseif comp==5     % Solución analítica de z2(xx,yy,tt)
     sol=exp(tt)*pi*(1 + xx)^2*sin(pi*xx)*sin(pi*yy) - exp(tt)*pi*cos(pi*xx)*cos(pi*yy)*sin(xx*yy);
 else               % Solución analítica de gamma(xx,yy,tt) (1st row & 2nd column)
     sol=(exp(tt)*(4*xx^3*yy^3 + (-1 + xx)*cos(1 - yy)*cos((1 - xx)*(1 - yy)) + sin(1 - yy)*sin((1 - xx)*(1 - yy))) -...
         exp(tt)*(-4*(1 - xx)^3*(1 - yy)^3 + cos(xx)*cos(xx*yy) - yy*sin(xx)*sin(xx*yy)))/2.;
 end

return
end
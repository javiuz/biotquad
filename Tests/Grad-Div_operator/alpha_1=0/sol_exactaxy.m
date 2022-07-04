function sol=sol_exactaxy(xx,yy,tt)

%       u(xx,yy,tt) = {e^t (x^2 + x^3 y^4 + Cos[1 - y] Sin[(1 - x) (1 - y)]), 
%                     e^t ((1 - y)^2 + (1 - x)^4 (1 - y)^3 + Cos[x y] Sin[x])}

%         u1(xx,yy) = e^t (x^2 + x^3 y^4 + Cos[1 - y] Sin[(1 - x) (1 - y)]);
%         u2(xx,yy) = e^t ((1 - y)^2 + (1 - x)^4 (1 - y)^3 + Cos[x y] Sin[x]);

 sol=[exp(tt)*(xx^2 + xx^3*yy^4 + cos(1-yy)*sin((1-xx)*(1-yy))),...
      exp(tt)*((1-yy)^2 + (1-xx)^4*(1-yy)^3 + cos(xx*yy)*sin(xx))]';

return
end
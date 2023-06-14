function sol=sol_exactaxy(xx,yy,tt)

%       u(xx,yy,tt) = {u1(xx,yy,tt, u2(xx,yy,tt}}

 sol=exp(tt)*[xx^2 + xx^3*yy^4 + cos(1 - yy)*sin((1 - xx)*(1 - yy)),...
             (1 - yy)^2 + (1 - xx)^4*(1 - yy)^3 + cos(xx*yy)*sin(xx)]';

return
end
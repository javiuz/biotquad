function sol=sol_exactaxy(xx,yy,~)

%       u(xx,yy,tt) = {u1(xx,yy,tt, u2(xx,yy,tt}}

 sol=[(1 - xx)*xx*yy,xx*(1 - yy)*yy]';

return
end
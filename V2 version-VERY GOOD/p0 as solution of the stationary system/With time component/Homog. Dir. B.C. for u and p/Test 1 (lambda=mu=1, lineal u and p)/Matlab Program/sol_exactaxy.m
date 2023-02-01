function sol=sol_exactaxy(xx,yy,tt)

%       u(xx,yy,tt) = {u1(xx,yy,tt, u2(xx,yy,tt}}

 sol=[tt*(1 - xx)*xx*(1 - yy)*yy,tt*(1 - xx)*xx*(1 - yy)*yy]';

return
end
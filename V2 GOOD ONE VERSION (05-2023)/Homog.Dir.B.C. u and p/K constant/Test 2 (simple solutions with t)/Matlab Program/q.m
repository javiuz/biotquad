function qq=q(xx,yy,tt)

global c0 perm

%       qq(xx,yy,tt) = q
    
qq=c0*(-1 + xx)*xx*(-1 + yy)*yy - 2*perm*tt*(-xx + xx^2 + (-1 + yy)*yy) + ...
    alpha*(-((-1 + yy)*yy) + xx^2*(-1 + 2*yy) + xx*(1 - 4*yy + 2*yy^2)); 
return
end
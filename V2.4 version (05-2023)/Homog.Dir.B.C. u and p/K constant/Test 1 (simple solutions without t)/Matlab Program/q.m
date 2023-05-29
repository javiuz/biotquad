function qq=q(xx,yy,tt)

global c0 perm

%       qq(xx,yy,tt) = q
    
qq=-2*perm*(-xx + xx^2 + (-1 + yy)*yy); 
return
end
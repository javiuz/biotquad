function qq=q0(xx,yy,tt)

global perm

%       qq(xx,yy,tt) = q
    
qq=-2*perm*(-xx + xx^2 + (-1 + yy)*yy); 
return
end
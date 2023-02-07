function qq_init=q0(xx,yy,tt)

global c0 perm

%       qq(xx,yy,tt) = q
    
qq_init=-2*perm*(-xx + xx^2 + (-1 + yy)*yy); 
return
end
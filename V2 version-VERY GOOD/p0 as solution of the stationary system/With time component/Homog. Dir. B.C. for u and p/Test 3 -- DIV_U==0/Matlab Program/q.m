function qq=q(xx,yy,tt)

global c0 perm alpha

%       qq(xx,yy,tt) = q
    
qq=c0*(-1 + xx)*xx*(-1 + yy)*yy - 2*perm*tt*(-xx + xx^2 + (-1 + yy)*yy); 
return
end
function qq=q(xx,yy,tt)

global c0 perm alpha

%       qq(xx,yy,tt) = q
    
qq=-2*perm*(-1 + xx)*xx - 2*perm*(-1 + yy)*yy + ...
    alpha*exp(tt)*(xx*(2 + 3*xx*yy^4) + (-1 + yy)*cos(1 - yy)*cos((-1 + xx)*(-1 + yy))) + ...
    c0*exp(tt)*(10 + cos(pi*yy)*sin(pi*xx)) + alpha*exp(tt)*...
    (2*(-1 + yy) - 3*(-1 + xx)^4*(-1 + yy)^2 - xx*sin(xx)*sin(xx*yy)); 
return
end
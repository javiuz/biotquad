function qq=q(xx,yy,tt)

global c0 perm alpha

%       qq(xx,yy,tt) = q
    
qq= exp(tt)*(2*alpha*(-1 + yy) - 3*alpha*(-1 + xx)^4*(-1 + yy)^2 + alpha*xx*(2 + 3*xx*yy^4) + ...
    alpha*(-1 + yy)*cos(1 - yy)*cos((-1 + xx)*(-1 + yy)) + 2*perm*pi^2*cos(pi*yy)*sin(pi*xx) + ...
    c0*(10 + cos(pi*yy)*sin(pi*xx)) - alpha*xx*sin(xx)*sin(xx*yy)); 
return
end
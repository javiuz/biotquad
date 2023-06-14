function qq=q(xx,yy,tt)

global c0 perm alpha

%       qq(xx,yy,tt) = q
    
qq= exp(tt)*(2*perm*pi^2*cos(pi*yy)*sin(pi*xx) + c0*(10 + cos(pi*yy)*sin(pi*xx))); 
return
end
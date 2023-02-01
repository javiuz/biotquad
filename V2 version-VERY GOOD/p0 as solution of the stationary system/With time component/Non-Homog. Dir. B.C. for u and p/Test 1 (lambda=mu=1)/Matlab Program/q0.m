function qq_init=q0(xx,yy,tt)

global c0 perm

%       qq(xx,yy,tt) = q
    
qq_init=2*exp(tt)*perm*pi^2*cos(pi*yy)*sin(pi*xx); 
return
end
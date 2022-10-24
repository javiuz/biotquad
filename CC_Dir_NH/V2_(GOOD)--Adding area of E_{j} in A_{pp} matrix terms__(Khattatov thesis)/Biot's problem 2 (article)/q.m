function qq=q(xx,yy,tt)

global c0 alpha

%       qq(xx,yy,tt) = q
    
qq=exp(tt)*(c0*(10 + cos(pi*yy)*sin(pi*xx)) +...
          alpha*(2*xx + 2*(-1+yy) - 3*(-1+xx)^4*(-1+yy)^2 + 3*xx^2*yy^4 +...
          (-1+yy)*cos(1-yy)*cos((-1+xx)*(-1+yy)) - xx*sin(xx)*sin(xx*yy)) +...
          pi*(sin(pi*xx)*(pi*(2+4*xx+2*xx^2+yy^2)*cos(pi*yy) + yy*cos(xx*yy)*sin(pi*yy)) -...
          cos(pi*xx)*(cos(pi*yy)*(2+2*xx+xx*cos(xx*yy)) - 2*pi*sin(pi*yy)*sin(xx*yy))));
       
return
end
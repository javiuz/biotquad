function qq=q(xx,yy,tt)

global c0 %alpha

%       qq(xx,yy,tt) = q
    
qq=exp(tt)*(10*c0 + (c0 + pi^2*(2 + 4*xx + 2*xx^2 + yy^2))*cos(pi*yy)*sin(pi*xx) + ...
          pi*yy*cos(xx*yy)*sin(pi*xx)*sin(pi*yy) + ...
          pi*cos(pi*xx)*(-(cos(pi*yy)*(2 + 2*xx + xx*cos(xx*yy))) + 2*pi*sin(pi*yy)*sin(xx*yy)));
       
return
end
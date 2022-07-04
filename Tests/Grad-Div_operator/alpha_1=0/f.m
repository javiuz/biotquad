function ff=f(xx,yy,k,tt,nu)

% global alpha

%       ff(xx,yy,t) = {f1,f2}'
    
if k==1     % Primera componente del vector ff (f1)          
    ff=(exp(tt)*(-((5+ sin(5*pi*xx)*sin(5*pi*yy))*...
        (-12*(-1 + xx)^3*(-1 + yy)^2 + 12*xx^3*yy^2 - xx*yy*cos(xx*yy)*sin(xx) +...
        2*(-1 + xx)*cos((-1 + xx)*(-1 + yy))*sin(1 - yy) - cos(1 - yy)*sin((-1 + xx)*(-1 + yy)) -...
        (-1 + xx)^2*cos(1 - yy)*sin((-1 + xx)*(-1 + yy)) - xx*cos(xx)*sin(xx*yy) - sin(xx)*sin(xx*yy))) -...
        5*pi*cos(5*pi*yy)*sin(5*pi*xx)*(-4*(-1 + xx)^3*(-1 + yy)^3 + 4*xx^3*yy^3 + ...
        (-1 + xx)*cos(1 - yy)*cos((-1 + xx)*(-1 + yy)) + cos(xx)*cos(xx*yy) + ...
        sin(1 - yy)*sin((-1 + xx)*(-1 + yy)) - yy*sin(xx)*sin(xx*yy)) - ...
        2*(5 + sin(5*pi*xx)*sin(5*pi*yy))*(2 + 6*xx*yy^4 - (-1 + yy)^2*cos(1 - yy)*sin((-1 + xx)*(-1 + yy)) + ...
        (nu*(2 - 12*(-1 + xx)^3*(-1 + yy)^2 + 6*xx*yy^4 - xx*yy*cos(xx*yy)*sin(xx) - ...
        (-1 + yy)^2*cos(1 - yy)*sin((-1 + xx)*(-1 + yy)) - xx*cos(xx)*sin(xx*yy) - sin(xx)*sin(xx*yy)))/...
        (1 - 2*nu)) - 10*pi*cos(5*pi*xx)*sin(5*pi*yy)*...
        (xx*(2+ 3*xx*yy^4) + (-1 + yy)*cos(1 - yy)*cos((-1 + xx)*(-1 + yy)) + ...
        (nu*(2*(-1 + yy) - 3*(-1 + xx)^4*(-1 + yy)^2 + xx*(2 + 3*xx*yy^4) + ...
        (-1 + yy)*cos(1 - yy)*cos((-1 + xx)*(-1 + yy)) - xx*sin(xx)*sin(xx*yy)))/(1 - 2*nu))))/(2.*(1+nu));
else        % Segunda componente del vector ff (f2)   
    ff=(exp(tt)*(-2*(2 - 6*(-1 + xx)^4*(-1 + yy) - xx^2*cos(xx*yy)*sin(xx) + ...
       (nu*(2 - 6*(-1 + xx)^4*(-1 + yy) + 12*xx^2*yy^3 + cos(1 - yy)*cos((-1 + xx)*(-1 + yy)) - ...
       xx^2*cos(xx*yy)*sin(xx) + (-1 + yy)*cos((-1 + xx)*(-1 + yy))*sin(1 - yy) - ...
       (-1 + xx)*(-1 + yy)*cos(1 - yy)*sin((-1 + xx)*(-1 + yy))))/(1 - 2*nu))*(5 + sin(5*pi*xx)*sin(5*pi*yy))...
       - (5 + sin(5*pi*xx)*sin(5*pi*yy))*(-12*(-1 + xx)^2*(-1 + yy)^3 + 12*xx^2*yy^3 + ...
       cos(1 - yy)*cos((-1 + xx)*(-1 + yy)) - cos(xx*yy)*sin(xx) - yy^2*cos(xx*yy)*sin(xx) +  ...
       (-1 + yy)*cos((-1 + xx)*(-1 + yy))*sin(1 - yy) -  ...
       (-1 + xx)*(-1 + yy)*cos(1 - yy)*sin((-1 + xx)*(-1 + yy)) - 2*yy*cos(xx)*sin(xx*yy)) -  ...
       5*pi*cos(5*pi*xx)*sin(5*pi*yy)*(-4*(-1 + xx)^3*(-1 + yy)^3 + 4*xx^3*yy^3 + ...
       (-1 + xx)*cos(1 - yy)*cos((-1 + xx)*(-1 + yy)) + cos(xx)*cos(xx*yy) +  ...
       sin(1 - yy)*sin((-1 + xx)*(-1 + yy)) - yy*sin(xx)*sin(xx*yy)) - ...
       10*pi*cos(5*pi*yy)*sin(5*pi*xx)*(2*(-1 + yy) - 3*(-1 + xx)^4*(-1 + yy)^2 - xx*sin(xx)*sin(xx*yy) + ...
       (nu*(2*(-1 + yy) - 3*(-1 + xx)^4*(-1 + yy)^2 + xx*(2 + 3*xx*yy^4) + ...
       (-1 + yy)*cos(1 - yy)*cos((-1 + xx)*(-1 + yy)) - xx*sin(xx)*sin(xx*yy)))/(1 - 2*nu))))/(2.*(1 + nu));
end
       
return
end
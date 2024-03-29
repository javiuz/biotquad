function ff=f(xx,yy,k,tt)

global alpha lambda mu

%       ff(xx,yy,t) = {f1,f2}'
    
if k==1     % Primera componente del vector ff (f1)          
%    ff=alpha*tt;
%    ff=1;
    ff=lambda*(-1 + xx*(2 - 4*yy) + 4*yy - 2*yy^2) - mu*(1 + 2*xx^2 + 4*xx*(-1 + yy) - 6*yy + 4*yy^2);
else        % Segunda componente del vector ff (f2)   
%    ff=0;
     ff=lambda*(-1 - 2*xx^2 - 4*xx*(-1 + yy) + 2*yy) - mu*(1 + 4*xx^2 - 4*yy + 2*yy^2 + xx*(-6 + 4*yy));
end
       
return
end
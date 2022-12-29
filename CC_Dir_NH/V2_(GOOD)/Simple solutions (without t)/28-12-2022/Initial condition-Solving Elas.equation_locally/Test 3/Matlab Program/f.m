function ff=f(xx,yy,k,tt)

global alpha lambda mu

%       ff(xx,yy,t) = {f1,f2}'
    
if k==1     % Primera componente del vector ff (f1)          
    ff=lambda*(-1 + 4*yy) + mu*(-1 + 6*yy) + alpha*yy*(1 - 2*xx*yy);
else        % Segunda componente del vector ff (f2)   
    ff=lambda*(-1 + 4*xx) + mu*(-1 + 6*xx) + alpha*xx*(1 - 2*xx*yy);
end
       
return
end
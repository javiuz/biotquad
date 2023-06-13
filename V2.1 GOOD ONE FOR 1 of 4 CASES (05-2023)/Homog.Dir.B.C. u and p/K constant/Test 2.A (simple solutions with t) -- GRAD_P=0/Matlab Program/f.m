function ff=f(xx,yy,k,tt)

global alpha lambda mu

%       ff(xx,yy,t) = {f1,f2}'
    
if k==1     % Primera componente del vector ff (f1)          
    ff= -(tt*(mu*(1 + 2*xx^2 + 4*xx*(-1 + yy) - 6*yy + 4*yy^2) + ...
        lambda*(1 - 4*yy + 2*yy^2 + xx*(-2 + 4*yy))));
else        % Segunda componente del vector ff (f2)   
    ff=-(tt*(lambda*(1 + 2*xx^2 + 4*xx*(-1 + yy) - 2*yy) +  ...
        mu*(1 - 6*xx + 4*xx^2 - 4*yy + 4*xx*yy + 2*yy^2)));
end
       
return
end
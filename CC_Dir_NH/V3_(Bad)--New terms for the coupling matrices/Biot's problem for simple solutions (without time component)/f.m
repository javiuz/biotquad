function ff=f(xx,yy,k,~)

global alpha

%       ff(xx,yy) = {f1,f2}'
    
if k==1     % Primera componente del vector ff (f1)          
    ff=alpha*yy;
else        % Segunda componente del vector ff (f2)   
    ff=alpha*xx;
end
       
return
end
function ff=f(xx,yy,k,tt)

global alpha

%       ff(xx,yy,t) = {f1,f2}'
    
if k==1     % Primera componente del vector ff (f1)          
%    ff=alpha*tt;
%    ff=1;
    ff=(-1 + 2*xx) *(-1 + yy)*yy;
else        % Segunda componente del vector ff (f2)   
%    ff=0;
    ff=(-1 + xx)* xx* (-1 + 2* yy);
end
       
return
end
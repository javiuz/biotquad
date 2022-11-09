function ff=f(~,~,k,tt)

global alpha

%       ff(xx,yy,t) = {f1,f2}'
    
if k==1     % Primera componente del vector ff (f1)          
%    ff=alpha*tt;
    ff=1;
else        % Segunda componente del vector ff (f2)   
    ff=0;
end
       
return
end
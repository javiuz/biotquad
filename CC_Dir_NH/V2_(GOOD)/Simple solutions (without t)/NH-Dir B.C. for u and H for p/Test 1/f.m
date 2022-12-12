function ff=f(xx,yy,k,tt)

global alpha lambda mu

%       ff(xx,yy,t) = {f1,f2}'
    
if k==1     % Primera componente del vector ff (f1)          
%    ff=alpha*tt;
%    ff=1;
    ff=-2*(lambda + 2*mu);
else        % Segunda componente del vector ff (f2)   
%    ff=0;
     ff=-2*(lambda + 2*mu);
end
       
return
end
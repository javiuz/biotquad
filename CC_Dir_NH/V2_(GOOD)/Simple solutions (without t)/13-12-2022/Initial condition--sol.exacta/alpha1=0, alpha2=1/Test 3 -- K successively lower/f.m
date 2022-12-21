function ff=f(xx,yy,k,tt)

global alpha lambda mu

%       ff(xx,yy,t) = {f1,f2}'
    
if k==1     % Primera componente del vector ff (f1)          
    ff=alpha;
%     ff=0;
else        % Segunda componente del vector ff (f2)   
    ff=0;
end
       
return
end
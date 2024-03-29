function indep1=build_indep_f_integral2(tt)

global NN alpha lambda mu x y

N=NN;

indep1=zeros(2*N*N,1);

% 1st component of the source term f(f1)   

fun1 = @(x,y) lambda + alpha*y -2*alpha*x.*y -alpha*y.^2 + 2*alpha*x.*y.^2 -...
              2*lambda*y.^2 -2*lambda*x + 4*lambda*x.*y + mu - 2*mu*x.^2 +...
              2*mu*y + 4*mu*x.*y - 4*mu*y.^2;
          
% 2nd component of the source term f(f2)   
          
fun2 = @(x,y) -lambda + 2*lambda*x.^2 + 2*lambda*y - 4*lambda*x.*y -mu + ...
              2*mu*x - 4*mu*x.*y + 2*mu*y.^2 + alpha*x - 4*mu*x - ...
              alpha*x.^2 + 4*mu*x.^2 - 2*alpha*x.*y + 2*alpha*y.*x.^2;


for j=1:N
    for i=1:N
        ind2=(i+(j-1)*N)*2;
        ind1=ind2-1;
        x1=x(i,j);
        y1=y(i,j);
        x2=x(i+1,j);
        y2=y(i+1,j);
        x3=x(i+1,j+1);
        y3=y(i+1,j+1);
        x4=x(i,j+1);
        y4=y(i,j+1);
        xmin=min(x1,x4);
        xmax=max(x2,x3);
        ymin=min(y1,y2);
        ymax=max(y3,y4);
        indep1(ind1)=integral2(fun1,xmin,xmax,ymin,ymax);
        indep1(ind2)=integral2(fun2,xmin,xmax,ymin,ymax);
%         indep1(ind1)=integral2(fun1,xmin,xmax,ymin,ymax,...
%                                'AbsTol',1e-12,'RelTol',1e-10);
%         indep1(ind2)=integral2(fun2,xmin,xmax,ymin,ymax,...
%                                'AbsTol',1e-12,'RelTol',1e-10);
    end
end
return
end
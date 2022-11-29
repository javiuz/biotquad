function indep2=build_indep_q_integral2(tt)

global NN alpha lambda mu perm x y

N=NN;

indep2=zeros(N*N,1);

% Source term q   

fun = @(x,y) 2*perm*x - 2*perm*x.^2 + 2*perm*y - 2*perm*y.^2;

for j=1:N
    for i=1:N
        ind=i+(j-1)*N;  
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
        indep2(ind)=integral2(fun,xmin,xmax,ymin,ymax);
%         indep2(ind)=integral2(fun,xmin,xmax,ymin,ymax,...
%                               'AbsTol',1e-12,'RelTol',1e-10);
    end
end
return
end
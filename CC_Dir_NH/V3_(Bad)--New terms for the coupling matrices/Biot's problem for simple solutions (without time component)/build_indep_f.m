function indep1=build_indep_f

global NN 

N=NN;

indep1=zeros(2*N*N,1);

for j=1:N
    for i=1:N
        ind2=(i+(j-1)*N)*2;
        ind1=ind2-1;        
        indep1(ind1)=int_f_simpson(i,j,1);
        indep1(ind2)=int_f_simpson(i,j,2);
    end
end
return
end

function intt=int_f_simpson(i,j,k)

intt=1/36*(g(i,j,0,0,k)+4*g(i,j,0.5,0,k)+g(i,j,1,0,k)+...
            4*(g(i,j,0,0.5,k)+4*g(i,j,0.5,0.5,k)+g(i,j,1,0.5,k))+...
            g(i,j,0,1,k)+4*g(i,j,0.5,1,k)+g(i,j,1,1,k));
            
return 
end

function gg=g(i,j,xx,yy,k)

global x y 

t1=area_triangulo(x(i,j+1),y(i,j+1),x(i,j),y(i,j),x(i+1,j),y(i+1,j));
t2=area_triangulo(x(i,j),y(i,j),x(i+1,j),y(i+1,j),x(i+1,j+1),y(i+1,j+1));
t4=area_triangulo(x(i+1,j+1),y(i+1,j+1),x(i,j+1),y(i,j+1),x(i,j),y(i,j));

jacob=2*t1+2*(t2-t1)*xx+2*(t4-t1)*yy;
gg=jacob*f(x(i,j)+(x(i+1,j)-x(i,j))*xx+(x(i,j+1)-x(i,j))*yy+(x(i+1,j+1)-x(i,j+1)-x(i+1,j)+x(i,j))*xx*yy,...
           y(i,j)+(y(i+1,j)-y(i,j))*xx+(y(i,j+1)-y(i,j))*yy+(y(i+1,j+1)-y(i,j+1)-y(i+1,j)+y(i,j))*xx*yy,k);

return
end

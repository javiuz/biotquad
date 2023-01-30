function indep2=build_indep_q0(tt)

global NN 

N=NN;

indep2=zeros(N*N,1);

for j=1:N
    for i=1:N
        ind=i+(j-1)*N;  
%         xx=(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1))/4;
%         yy=(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1))/4;
%         E=sin(5*pi*xx)*sin(5*pi*yy)+5;
%         lambda=E*c_lambda;
%         mu=E*c_mu;      
        indep2(ind)=int_f_simpson(i,j,tt);
    end
end
return
end

function intt=int_f_simpson(i,j,tt)

intt=1/36*(g(i,j,0,0,tt)+4*g(i,j,0.5,0,tt)+g(i,j,1,0,tt)+...
            4*(g(i,j,0,0.5,tt)+4*g(i,j,0.5,0.5,tt)+g(i,j,1,0.5,tt))+...
            g(i,j,0,1,tt)+4*g(i,j,0.5,1,tt)+g(i,j,1,1,tt));
            
return 
end

function gg=g(i,j,xx,yy,tt)

global x y

t1=area_triangulo(x(i,j+1),y(i,j+1),x(i,j),y(i,j),x(i+1,j),y(i+1,j));
t2=area_triangulo(x(i,j),y(i,j),x(i+1,j),y(i+1,j),x(i+1,j+1),y(i+1,j+1));
t4=area_triangulo(x(i+1,j+1),y(i+1,j+1),x(i,j+1),y(i,j+1),x(i,j),y(i,j));

jacob=2*t1+2*(t2-t1)*xx+2*(t4-t1)*yy;
gg=jacob*q0(x(i,j)+(x(i+1,j)-x(i,j))*xx+(x(i,j+1)-x(i,j))*yy+(x(i+1,j+1)-x(i,j+1)-x(i+1,j)+x(i,j))*xx*yy,...
           y(i,j)+(y(i+1,j)-y(i,j))*xx+(y(i,j+1)-y(i,j))*yy+(y(i+1,j+1)-y(i,j+1)-y(i+1,j)+y(i,j))*xx*yy,tt);

return
end

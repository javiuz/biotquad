function indep1=build_indep_f_gauss_lobatto(tt)

global NN 

N=NN;

indep1=zeros(2*N*N,1);

for j=1:N
    for i=1:N
        ind2=(i+(j-1)*N)*2;
        ind1=ind2-1;    
%         xx=(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1))/4;
%         yy=(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1))/4;
%         E=sin(5*pi*xx)*sin(5*pi*yy)+5;
%         lambda=E*c_lambda;
%         mu=E*c_mu;      
        indep1(ind1)=int_f_gauss_lobatto(i,j,1,tt);
        indep1(ind2)=int_f_gauss_lobatto(i,j,2,tt);
    end
end
return
end

function intt=int_f_gauss_lobatto(i,j,k,tt)

x1=0;
x2=0.2763932;
x3=0.7236068;
x4=1;
w1=1/12;
w2=5/12;

intt=w1*(w1*g(i,j,x1,x1,k,tt)+w2*g(i,j,x2,x1,k,tt)+w2*g(i,j,x3,x1,k,tt)+w1*g(i,j,x4,x1,k,tt))+...
    w2*(w1*g(i,j,x1,x2,k,tt)+w2*g(i,j,x2,x2,k,tt)+w2*g(i,j,x3,x2,k,tt)+w1*g(i,j,x4,x2,k,tt))+...
    w2*(w1*g(i,j,x1,x3,k,tt)+w2*g(i,j,x2,x3,k,tt)+w2*g(i,j,x3,x3,k,tt)+w1*g(i,j,x4,x3,k,tt))+...
    w1*(w1*g(i,j,x1,x4,k,tt)+w2*g(i,j,x2,x4,k,tt)+w2*g(i,j,x3,x4,k,tt)+w1*g(i,j,x4,x4,k,tt));
            
return 
end

function gg=g(i,j,xx,yy,k,tt)

global x y 

t1=area_triangulo(x(i,j+1),y(i,j+1),x(i,j),y(i,j),x(i+1,j),y(i+1,j));
t2=area_triangulo(x(i,j),y(i,j),x(i+1,j),y(i+1,j),x(i+1,j+1),y(i+1,j+1));
t4=area_triangulo(x(i+1,j+1),y(i+1,j+1),x(i,j+1),y(i,j+1),x(i,j),y(i,j));

jacob=2*t1+2*(t2-t1)*xx+2*(t4-t1)*yy;
gg=jacob*f(x(i,j)+(x(i+1,j)-x(i,j))*xx+(x(i,j+1)-x(i,j))*yy+(x(i+1,j+1)-x(i,j+1)-x(i+1,j)+x(i,j))*xx*yy,...
           y(i,j)+(y(i+1,j)-y(i,j))*xx+(y(i,j+1)-y(i,j))*yy+(y(i+1,j+1)-y(i,j+1)-y(i+1,j)+y(i,j))*xx*yy,k,tt);

return
end
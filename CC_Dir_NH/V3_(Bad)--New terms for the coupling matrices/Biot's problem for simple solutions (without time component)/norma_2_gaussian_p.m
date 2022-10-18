function norm=norma_2_gaussian_p(pp,i,j)

x1=(5-sqrt(15))/10;
x2=1/2;
x3=(5+sqrt(15))/10;
w1=5/18;
w2=4/9;

norm=w1*(w1*g_gaussian(pp,i,j,x1,x1)+w2*g_gaussian(pp,i,j,x2,x1)+...
    w1*g_gaussian(pp,i,j,x3,x1))+w2*(w1*g_gaussian(pp,i,j,x1,x2)+...
    w2*g_gaussian(pp,i,j,x2,x2)+w1*g_gaussian(pp,i,j,x3,x2))+...
    w1*(w1*g_gaussian(pp,i,j,x1,x3)+w2*g_gaussian(pp,i,j,x2,x3)+...
    w1*g_gaussian(pp,i,j,x3,x3));

return
end 

function gg=g_gaussian(pp,i,j,xx,yy)

global x y

t1=area_triangulo(x(i,j+1),y(i,j+1),x(i,j),y(i,j),x(i+1,j),y(i+1,j));
t2=area_triangulo(x(i,j),y(i,j),x(i+1,j),y(i+1,j),x(i+1,j+1),y(i+1,j+1));
t4=area_triangulo(x(i+1,j+1),y(i+1,j+1),x(i,j+1),y(i,j+1),x(i,j),y(i,j));

jacob=2*t1+2*(t2-t1)*xx+2*(t4-t1)*yy;
gg=jacob*(sol_exactax(x(i,j)+(x(i+1,j)-x(i,j))*xx+(x(i,j+1)-x(i,j))*yy+...
        (x(i+1,j+1)-x(i,j+1)-x(i+1,j)+x(i,j))*xx*yy,...
        y(i,j)+(y(i+1,j)-y(i,j))*xx+(y(i,j+1)-y(i,j))*yy+...
        (y(i+1,j+1)-y(i,j+1)-y(i+1,j)+y(i,j))*xx*yy,3)-pp)^2;

return
end 

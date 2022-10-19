function norm=norma_2_gaussian_vel(i,j,vel,t)

% sigma is a vector with 16 components

x1=(5-sqrt(15))/10;
x2=1/2;
x3=(5+sqrt(15))/10;
w1=5/18;
w2=4/9;

norm=w1*(w1*g_gaussian_vel(i,j,vel,x1,x1,t)+w2*g_gaussian_vel(i,j,vel,x2,x1,t)+...
    w1*g_gaussian_vel(i,j,vel,x3,x1,t))+w2*(w1*g_gaussian_vel(i,j,vel,x1,x2,t)+...
    w2*g_gaussian_vel(i,j,vel,x2,x2,t)+w1*g_gaussian_vel(i,j,vel,x3,x2,t))+...
    w1*(w1*g_gaussian_vel(i,j,vel,x1,x3,t)+w2*g_gaussian_vel(i,j,vel,x2,x3,t)+...
    w1*g_gaussian_vel(i,j,vel,x3,x3,t));

return
end 

function gg=g_gaussian_vel(i,j,vel,xx,yy,tt)

global x y 

t1=area_triangulo(x(i,j+1),y(i,j+1),x(i,j),y(i,j),x(i+1,j),y(i+1,j));
t2=area_triangulo(x(i,j),y(i,j),x(i+1,j),y(i+1,j),x(i+1,j+1),y(i+1,j+1));
t4=area_triangulo(x(i+1,j+1),y(i+1,j+1),x(i,j+1),y(i,j+1),x(i,j),y(i,j));

jacob=2*t1+2*(t2-t1)*xx+2*(t4-t1)*yy;

xg=x(i,j)+(x(i+1,j)-x(i,j))*xx+(x(i,j+1)-x(i,j))*yy+...
        (x(i+1,j+1)-x(i,j+1)-x(i+1,j)+x(i,j))*xx*yy;
yg=y(i,j)+(y(i+1,j)-y(i,j))*xx+(y(i,j+1)-y(i,j))*yy+...
        (y(i+1,j+1)-y(i,j+1)-y(i+1,j)+y(i,j))*xx*yy;
    
df=[x(i+1,j)-x(i,j) x(i,j+1)-x(i,j);y(i+1,j)-y(i,j) y(i,j+1)-y(i,j)]+...
[(x(i+1,j+1)-x(i,j+1)-x(i+1,j)+x(i,j))*yy (x(i+1,j+1)-x(i,j+1)-x(i+1,j)+x(i,j))*xx;...
 (y(i+1,j+1)-y(i,j+1)-y(i+1,j)+y(i,j))*yy (y(i+1,j+1)-y(i,j+1)-y(i+1,j)+y(i,j))*xx];

[a1,b1,g1,a2,b2,g2,r,s]=compute_vel_coef(vel);
sol_nca_1=a1*xx+b1*yy+g1+r*xx^2+2*s*xx*yy;
sol_nca_2=a2*xx+b2*yy+g2-2*r*xx*yy-s*yy^2;

zg=[sol_exactax(xg,yg,tt,4);sol_exactax(xg,yg,tt,5)];
zhg=[sol_nca_1;sol_nca_2];

gg=jacob*(zg-(1/jacob)*df*zhg)'*(zg-(1/jacob)*df*zhg);

return
end 

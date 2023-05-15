function norm=norma_2_gauss_lobatto_rot(g1,g2,g3,g4,i,j,t)

x1=0;
x2=0.2763932;
x3=0.7236068;
x4=1;
w1=1/12;
w2=5/12;

norm=w1*(w1*g_gaussian_rot(g1,g2,g3,g4,i,j,x1,x1,t)+w2*g_gaussian_rot(g1,g2,g3,g4,i,j,x2,x1,t)+...
    w2*g_gaussian_rot(g1,g2,g3,g4,i,j,x3,x1,t)+...
    w1*g_gaussian_rot(g1,g2,g3,g4,i,j,x4,x1,t))+w2*(w1*g_gaussian_rot(g1,g2,g3,g4,i,j,x1,x2,t)+...
    w2*g_gaussian_rot(g1,g2,g3,g4,i,j,x2,x2,t)+w2*g_gaussian_rot(g1,g2,g3,g4,i,j,x3,x2,t)+w1*g_gaussian_rot(g1,g2,g3,g4,i,j,x4,x2,t))+...
    +w2*(w1*g_gaussian_rot(g1,g2,g3,g4,i,j,x1,x3,t)+...
    w2*g_gaussian_rot(g1,g2,g3,g4,i,j,x2,x3,t)+w2*g_gaussian_rot(g1,g2,g3,g4,i,j,x3,x3,t)+w1*g_gaussian_rot(g1,g2,g3,g4,i,j,x4,x3,t))+...
    w1*(w1*g_gaussian_rot(g1,g2,g3,g4,i,j,x1,x4,t)+w2*g_gaussian_rot(g1,g2,g3,g4,i,j,x2,x4,t)+...
    w2*g_gaussian_rot(g1,g2,g3,g4,i,j,x3,x4,t)+w1*g_gaussian_rot(g1,g2,g3,g4,i,j,x4,x4,t));

return
end 

function gg=g_gaussian_rot(g1,g2,g3,g4,i,j,xx,yy,tt)

global x y

t1=area_triangulo(x(i,j+1),y(i,j+1),x(i,j),y(i,j),x(i+1,j),y(i+1,j));
t2=area_triangulo(x(i,j),y(i,j),x(i+1,j),y(i+1,j),x(i+1,j+1),y(i+1,j+1));
t4=area_triangulo(x(i+1,j+1),y(i+1,j+1),x(i,j+1),y(i,j+1),x(i,j),y(i,j));

sol_nca=g1+(g2-g1)*xx+(g4-g1)*yy+(g1-g2+g3-g4)*xx*yy;

jacob=2*t1+2*(t2-t1)*xx+2*(t4-t1)*yy;
gg=jacob*((sol_exactax(x(i,j)+(x(i+1,j)-x(i,j))*xx+(x(i,j+1)-x(i,j))*yy+...
        (x(i+1,j+1)-x(i,j+1)-x(i+1,j)+x(i,j))*xx*yy,...
        y(i,j)+(y(i+1,j)-y(i,j))*xx+(y(i,j+1)-y(i,j))*yy+...
        (y(i+1,j+1)-y(i,j+1)-y(i+1,j)+y(i,j))*xx*yy,tt,6)-sol_nca)^2);

return
end 

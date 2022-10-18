function norm=norma_2_gaussian_stress(i,j,sigma,nu)

% sigma is a vector with 16 components

x1=(5-sqrt(15))/10;
x2=1/2;
x3=(5+sqrt(15))/10;
w1=5/18;
w2=4/9;

norm=w1*(w1*g_gaussian_stress(i,j,sigma,x1,x1,nu)+w2*g_gaussian_stress(i,j,sigma,x2,x1,nu)+...
    w1*g_gaussian_stress(i,j,sigma,x3,x1,nu))+w2*(w1*g_gaussian_stress(i,j,sigma,x1,x2,nu)+...
    w2*g_gaussian_stress(i,j,sigma,x2,x2,nu)+w1*g_gaussian_stress(i,j,sigma,x3,x2,nu))+...
    w1*(w1*g_gaussian_stress(i,j,sigma,x1,x3,nu)+w2*g_gaussian_stress(i,j,sigma,x2,x3,nu)+...
    w1*g_gaussian_stress(i,j,sigma,x3,x3,nu));

return
end 

function gg=g_gaussian_stress(i,j,sigma,xx,yy,nu)

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

[a1,b1,g1,a2,b2,g2,r,s]=compute_stress_row(sigma,1);
sol_nca_11=a1*xx+b1*yy+g1+r*xx^2+2*s*xx*yy;
sol_nca_12=a2*xx+b2*yy+g2-2*r*xx*yy-s*yy^2;

[a1,b1,g1,a2,b2,g2,r,s]=compute_stress_row(sigma,2);
sol_nca_21=a1*xx+b1*yy+g1+r*xx^2+2*s*xx*yy;
sol_nca_22=a2*xx+b2*yy+g2-2*r*xx*yy-s*yy^2;

sg=[sol_exactax_sigma(xg,yg,nu,1,1) sol_exactax_sigma(xg,yg,nu,1,2);sol_exactax_sigma(xg,yg,nu,2,1) sol_exactax_sigma(xg,yg,nu,2,2)];
shg=[sol_nca_11 sol_nca_12;sol_nca_21 sol_nca_22];


gg=jacob*scalar_prod_matrices(sg-1/jacob*shg*(df'),sg-1/jacob*shg*(df'));

return
end 

function s=scalar_prod_matrices(m,n)

s=m(1,1)*n(1,1)+m(1,2)*n(1,2)+m(2,1)*n(2,1)+m(2,2)*n(2,2);
end

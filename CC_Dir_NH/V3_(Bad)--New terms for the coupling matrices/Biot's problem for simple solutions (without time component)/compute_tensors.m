function [sigma,sx1,sx2,sy1,sy2]=compute_tensors(u,paux,gaux)

global NN x y alpha lambda mu

N=NN;

sigma=zeros(8*N*(N+1),1);

p=reshape(paux,N,N);
gamma=reshape(gaux,N+1,N+1);

sx1=zeros(N+1,N+1,2);
sx2=sx1;
sy1=sx1;
sy2=sx1;

coef_denom=16*mu*(lambda + mu);
coef_d=alpha/coef_denom;

% South-West corner node
i=1;
j=1;
ind2u=(i+(j-1)*N)*2;
ind1u=ind2u-1;
vdim=4;
pdim=1;
pp=zeros(pdim,1);
ld=1;

x1=x(i,j);  
y1=y(i,j);            
x2=x(i+1,j);
y2=y(i+1,j);
% x3=x(i+1,j+1);
% y3=y(i+1,j+1);
x4=x(i,j+1);
y4=y(i,j+1);

Jer1=abs(x4*(y1 - y2) + x1*(y2 - y4) + x2*(-y1 + y4));
denom_elas=coef_denom*Jer1;
denom_d=coef_d/Jer1;

% Local matrix (A sigma, sigma)
a=zeros(vdim,vdim);
a(1,1)=(lambda + 2*mu)*(x1 - x4)^2 + 2*(lambda + mu)*(y1 - y4)^2;
a(1,2)=lambda*(-x1 + x4)*(y1 - y4);
a(1,3)=(lambda + 2*mu)*(x1 - x2)*(x1 - x4) + 2*(lambda + mu)*(y1 - y2)*(y1 - y4);
a(1,4)=lambda*(-x1 + x4)*(y1 - y2);
a(2,1)=a(1,2);
a(2,2)=2*(lambda + mu)*(x1 - x4)^2 + (lambda + 2*mu)*(y1 - y4)^2;
a(2,3)=lambda*(-x1 + x2)*(y1 - y4);
a(2,4)=2*(lambda + mu)*(x1 - x2)*(x1 - x4) + (lambda + 2*mu)*(y1 - y2)*(y1 - y4);
a(3,1)=a(1,3);
a(3,2)=a(2,3);
a(3,3)=(lambda + 2*mu)*(x1 - x2)^2 + 2*(lambda + mu)*(y1 - y2)^2;
a(3,4)=lambda*(-x1 + x2)*(y1 - y2);
a(4,1)=a(1,4);
a(4,2)=a(2,4);
a(4,3)=a(3,4);
a(4,4)=2*(lambda + mu)*(x1 - x2)^2 + (lambda + 2*mu)*(y1 - y2)^2;
a=(1/denom_elas)*a;

% Local matrix (A sigma, u)^t
b=1/2*[-1 0;0 -1;-1 0;0 -1];

% Local matrix (A sigma, gamma)^t
c=1/4*[y4-y1 x1-x4 y2-y1 x1-x2]';

% Local matrix (A sigma, p)^t
d=zeros(vdim,1);
d(1,1)=(lambda + 2*mu)*(x1 - x2)*(x1 - x4) - lambda*(x1 - x4)*(y1 - y4) + 2*(lambda + mu)*(y1 - y2)*(y1 - y4);
d(2,1)=2*(lambda + mu)*(x1 - x4)^2 + lambda*(-x1 + x2)*(y1 - y4) + (lambda + 2*mu)*(y1 - y4)^2;
d(3,1)=(lambda + 2*mu)*(x1 - x2)^2 + 2*(lambda + mu)*(y1 - y2)^2 + lambda*(-x1 + x2)*(y1 - y4);
d(4,1)=2*(lambda + mu)*(x1 - x2)*(x1 - x4) + lambda*(-x1 + x2)*(y1 - y2) + (lambda + 2*mu)*(y1 - y2)*(y1 - y4);
d=denom_d*d;

% CC.D en el nodo S-W, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
    % frontera Sur
        Pg1u_L1=int_dir_simpson(i,j,i+1,j,1);
        Pg2u_L1=int_dir_simpson(i,j,i+1,j,2);
    % frontera Oeste
        Pg1u_L4=int_dir_simpson(i,j,i,j+1,1);
        Pg2u_L4=int_dir_simpson(i,j,i,j+1,2);      
Pgu=(1/2)*[-Pg1u_L1;-Pg2u_L1;-Pg1u_L4;-Pg2u_L4]; 

pp(1,1)=p(i,j);
gg=gamma(i,j);

ss=a\(Pgu-b*u(ind1u:ind2u)-c*gg-d*pp);
sigma(ld:vdim)=ss;
ld=vdim;

sy1(i,j,1)=ss(1,1);
sy2(i,j,1)=ss(2,1);
sx1(i,j,2)=ss(3,1);
sx2(i,j,2)=ss(4,1);

% South nodes (j=1)
vdim=6;
pdim=2;
pp=zeros(pdim,1);

% Matrix b is the same for every South node
b=(1/2)*[0 0 -1 0;0 0 0 -1;1 0 -1 0;0 1 0 -1;-1 0 0 0;0 -1 0 0];

a=zeros(vdim,vdim);
d=zeros(vdim,2);

for i=2:N

    ind2u=(i-1+(j-1)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i+(j-1)*N)*2;
%     ind3u=ind4u-1;

    x1=x(i-1,j);
    y1=y(i-1,j);
    x2=x(i,j);
    y2=y(i,j);
    x3=x(i,j+1);
    y3=y(i,j+1);
%     x4=x(i-1,j+1);
%     y4=y(i-1,j+1);
    x5=x(i+1,j);
    y5=y(i+1,j);
%     x6=x(i+1,j+1);
%     y6=y(i+1,j+1);
    
    JE1r2=abs(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3));
    denom_elas1=coef_denom*JE1r2;
    denom_d1=coef_d/JE1r2;
    JE2r1=abs(x5*(-y2 + y3) + x3*(y2 - y5) + x2*(-y3 + y5));
    denom_elas2=coef_denom*JE2r1;
    denom_d2=coef_d/JE2r1;
    
    a(1,1)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/denom_elas2;
    a(1,2)=-((lambda*(x2 - x3)*(y2 - y3))/denom_elas2);
    a(1,3)=((lambda + 2*mu)*(x2 - x3)*(x2 - x5) + 2*(lambda + mu)*(y2 - y3)*(y2 - y5))/denom_elas2;
    a(1,4)=-((lambda*(x2 - x3)*(y2 - y5))/denom_elas2);
%     a(1,5)=0;
%     a(1,6)=0;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/denom_elas2;
    a(2,3)=-((lambda*(x2 - x5)*(y2 - y3))/denom_elas2);
    a(2,4)=(2*(lambda + mu)*(x2 - x3)*(x2 - x5) + (lambda + 2*mu)*(y2 - y3)*(y2 - y5))/denom_elas2;
%     a(2,5)=0;
%     a(2,6)=0;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=((lambda + 2*mu)*(x1 - x2)^2 + 2*(lambda + mu)*(y1 - y2)^2)/denom_elas1 + ...
           ((lambda + 2*mu)*(x2 - x5)^2 + 2*(lambda + mu)*(y2 - y5)^2)/denom_elas2;
    a(3,4)=-((lambda*(x1 - x2)*(y1 - y2))/denom_elas1) - (lambda*(x2 - x5)*(y2 - y5))/denom_elas2;
    a(3,5)=((lambda + 2*mu)*(x1 - x2)*(x2 - x3) + 2*(lambda + mu)*(y1 - y2)*(y2 - y3))/denom_elas1;
    a(3,6)=-((lambda*(x1 - x2)*(y2 - y3))/denom_elas1);
    a(4,1)=a(1,4);
    a(4,2)=a(2,4);
    a(4,3)=a(3,4);
    a(4,4)=(2*(lambda + mu)*(x1 - x2)^2 + (lambda + 2*mu)*(y1 - y2)^2)/denom_elas1 + ...
           (2*(lambda + mu)*(x2 - x5)^2 + (lambda + 2*mu)*(y2 - y5)^2)/denom_elas2;
    a(4,5)=-((lambda*(x2 - x3)*(y1 - y2))/denom_elas1);
    a(4,6)=(2*(lambda + mu)*(x1 - x2)*(x2 - x3) + (lambda + 2*mu)*(y1 - y2)*(y2 - y3))/denom_elas1;
%     a(5,1)=0;
%     a(5,2)=0;
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/denom_elas1;
    a(5,6)=-((lambda*(x2 - x3)*(y2 - y3))/denom_elas1);
%     a(6,1)=0;
%     a(6,2)=0;
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/denom_elas1;
    
    c=(1/4)*[y3-y2 x2-x3 y5-y1 x1-x5 y3-y2 x2-x3]';
    
%     d(1,1)=0;
%     d(2,1)=0;
    d(3,1)=(lambda + 2*mu)*(x1 - x2)^2 + 2*(lambda + mu)*(y1 - y2)^2 + lambda*(-x1 + x2)*(y2 - y3);
    d(4,1)=2*(lambda + mu)*(x1 - x2)*(x2 - x3) + lambda*(-x1 + x2)*(y1 - y2) + (lambda + 2*mu)*(y1 - y2)*(y2 - y3);
    d(5,1)=(lambda + 2*mu)*(x1 - x2)*(x2 - x3) + lambda*(-x2 + x3)*(y2 - y3) + 2*(lambda + mu)*(y1 - y2)*(y2 - y3);
    d(6,1)=2*(lambda + mu)*(x2 - x3)^2 - lambda*(x1 - x2)*(y2 - y3) + (lambda + 2*mu)*(y2 - y3)^2;
    d(1,2)=(lambda + 2*mu)*(x2 - x3)*(x2 - x5) + lambda*(-x2 + x3)*(y2 - y3) + 2*(lambda + mu)*(y2 - y3)*(y2 - y5);
    d(2,2)=2*(lambda + mu)*(x2 - x3)^2 + lambda*(-x2 + x5)*(y2 - y3) + (lambda + 2*mu)*(y2 - y3)^2;
    d(3,2)=(lambda + 2*mu)*(x2 - x5)^2 + lambda*(-x2 + x5)*(y2 - y3) + 2*(lambda + mu)*(y2 - y5)^2;
    d(4,2)=2*(lambda + mu)*(x2 - x3)*(x2 - x5) - lambda*x2*(y2 - y5) + (lambda*x5 + (lambda + 2*mu)*(y2 - y3))*(y2 - y5);
%     d(5,2)=0;
%     d(6,2)=0;

    d(:,1)=denom_d1*d(:,1);
    d(:,2)=denom_d2*d(:,2);
    
    % CC.D en el nodo Sur, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
        % frontera Sur 1
            Pg1u_L2=int_dir_simpson(i,j,i+1,j,1);
            Pg2u_L2=int_dir_simpson(i,j,i+1,j,2);
        % frontera Sur 2
            Pg1u_L1=int_dir_simpson(i-1,j,i,j,1);
            Pg2u_L1=int_dir_simpson(i-1,j,i,j,2);
        Pgu=(1/2)*[-Pg1u_L2;-Pg2u_L2;0;0;-Pg1u_L1;-Pg2u_L1];
    
    pp(1,1)=p(i-1,j);
    pp(2,1)=p(i,j);
    gg=gamma(i,j);
    
    ss=a\(Pgu-b*u(ind1u:ind4u)-c*gg-d*pp);
    sigma(ld+1:ld+vdim)=ss;
    ld=ld+vdim;

    sy1(i,j,1)=ss(1,1);
    sy2(i,j,1)=ss(2,1);
    sx1(i,j,2)=ss(3,1);
    sx2(i,j,2)=ss(4,1);
    sy1(i,j,2)=ss(5,1);
    sy2(i,j,2)=ss(6,1); 
end 

% South-East corner node (j=1)
i=N+1;
ind2u=(i-1+(j-1)*N)*2;
ind1u=ind2u-1;
vdim=4;
pdim=1;
pp=zeros(pdim,1);

x1=x(i-1,j);  
y1=y(i-1,j);            
x2=x(i,j);
y2=y(i,j);
x3=x(i,j+1);
y3=y(i,j+1);
% x4=x(i-1,j+1);
% y4=y(i-1,j+1);

Jer2=abs(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3));
denom_elas=coef_denom*Jer2;
denom_d=coef_d/Jer2;

a=zeros(vdim,vdim);
a(1,1)=(lambda + 2*mu)*(x1 - x2)^2 + 2*(lambda + mu)*(y1 - y2)^2;
a(1,2)=lambda*(-x1 + x2)*(y1 - y2);
a(1,3)=(lambda + 2*mu)*(x1 - x2)*(x2 - x3) + 2*(lambda + mu)*(y1 - y2)*(y2 - y3);
a(1,4)=lambda*(-x1 + x2)*(y2 - y3);
a(2,1)=a(1,2);
a(2,2)=2*(lambda + mu)*(x1 - x2)^2 + (lambda + 2*mu)*(y1 - y2)^2;
a(2,3)=lambda*(-x2 + x3)*(y1 - y2);
a(2,4)=2*(lambda + mu)*(x1 - x2)*(x2 - x3) + (lambda + 2*mu)*(y1 - y2)*(y2 - y3);
a(3,1)=a(1,3);
a(3,2)=a(2,3);
a(3,3)=(lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2;
a(3,4)=lambda*(-x2 + x3)*(y2 - y3);
a(4,1)=a(1,4);
a(4,2)=a(2,4);
a(4,3)=a(3,4);
a(4,4)=2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2;
a=(1/denom_elas)*a;

b=(1/2)*[1 0;0 1;-1 0;0 -1];

c=(1/4)*[y2-y1 x1-x2 y3-y2 x2-x3]';

d=zeros(vdim,1);
d(1,1)=(lambda + 2*mu)*(x1 - x2)^2 + 2*(lambda + mu)*(y1 - y2)^2 + lambda*(-x1 + x2)*(y2 - y3);
d(2,1)=2*(lambda + mu)*(x1 - x2)*(x2 - x3) + lambda*(-x1 + x2)*(y1 - y2) + (lambda + 2*mu)*(y1 - y2)*(y2 - y3);
d(3,1)=(lambda + 2*mu)*(x1 - x2)*(x2 - x3) + lambda*(-x2 + x3)*(y2 - y3) + 2*(lambda + mu)*(y1 - y2)*(y2 - y3);
d(4,1)=2*(lambda + mu)*(x2 - x3)^2 - lambda*(x1 - x2)*(y2 - y3) + (lambda + 2*mu)*(y2 - y3)^2;
d=denom_d*d;

% CC.D en el nodo S-E, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
    % frontera Este
        Pg1u_L7=int_dir_simpson(i,j,i,j+1,1);
        Pg2u_L7=int_dir_simpson(i,j,i,j+1,2);
    % frontera Sur
        Pg1u_L3=int_dir_simpson(i-1,j,i,j,1);
        Pg2u_L3=int_dir_simpson(i-1,j,i,j,2);
    Pgu=(1/2)*[Pg1u_L7;Pg2u_L7;-Pg1u_L3;-Pg2u_L3];
    
pp(1,1)=p(i-1,j);
gg=gamma(i,j);

ss=a\(Pgu-b*u(ind1u:ind2u)-c*gg-d*pp);
sigma(ld+1:ld+vdim)=ss;
ld=ld+vdim;

sx1(i,j,2)=ss(1,1);
sx2(i,j,2)=ss(2,1);
sy1(i,j,2)=ss(3,1);
sy2(i,j,2)=ss(4,1);

for j=2:N

    % West nodes
    i=1;
    ind2u=(i+(j-2)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i+(j-1)*N)*2;
    ind3u=ind4u-1;
    vdim=6;
    pdim=2;
    pp=zeros(pdim,1);
    
    x1=x(i,j-1);
    y1=y(i,j-1);
%     x2=x(i+1,j-1);
%     y2=y(i+1,j-1);
    x3=x(i+1,j);
    y3=y(i+1,j);
    x4=x(i,j);
    y4=y(i,j);
%     x5=x(i+1,j+1);
%     y5=y(i+1,j+1);
    x6=x(i,j+1);
    y6=y(i,j+1);

    JE1r4=abs(x4*(y1 - y3) + x1*(y3 - y4) + x3*(-y1 + y4));
    denom_elas1=coef_denom*JE1r4;
    denom_d1=coef_d/JE1r4;
    JE2r1=abs(x6*(-y3 + y4) + x4*(y3 - y6) + x3*(-y4 + y6));
    denom_elas2=coef_denom*JE2r1;
    denom_d2=coef_d/JE2r1;
    
    a=zeros(vdim,vdim);    
    a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/denom_elas1;
    a(1,2)=-((lambda*(x3 - x4)*(y3 - y4))/denom_elas1);
    a(1,3)=-(((lambda + 2*mu)*(x1 - x4)*(x3 - x4) + 2*(lambda + mu)*(y1 - y4)*(y3 - y4))/denom_elas1);
    a(1,4)=(lambda*(x3 - x4)*(y1 - y4))/denom_elas1;
%     a(1,5)=0;
%     a(1,6)=0;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/denom_elas1;
    a(2,3)=(lambda*(x1 - x4)*(y3 - y4))/denom_elas1;
    a(2,4)=-((2*(lambda + mu)*(x1 - x4)*(x3 - x4) + (lambda + 2*mu)*(y1 - y4)*(y3 - y4))/denom_elas1);
%     a(2,5)=0;
%     a(2,6)=0;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=((lambda + 2*mu)*(x1 - x4)^2 + 2*(lambda + mu)*(y1 - y4)^2)/denom_elas1 + ...
           ((lambda + 2*mu)*(x4 - x6)^2 + 2*(lambda + mu)*(y4 - y6)^2)/denom_elas2;
    a(3,4)=-((lambda*(x1 - x4)*(y1 - y4))/denom_elas1) - (lambda*(x4 - x6)*(y4 - y6))/denom_elas2;
    a(3,5)=-(((lambda + 2*mu)*(x3 - x4)*(x4 - x6) + 2*(lambda + mu)*(y3 - y4)*(y4 - y6))/denom_elas2);
    a(3,6)=(lambda*(x4 - x6)*(y3 - y4))/denom_elas2;
    a(4,1)=a(1,4);
    a(4,2)=a(2,4);
    a(4,3)=a(3,4);
    a(4,4)=(2*(lambda + mu)*(x1 - x4)^2 + (lambda + 2*mu)*(y1 - y4)^2)/denom_elas1 + ...
           (2*(lambda + mu)*(x4 - x6)^2 + (lambda + 2*mu)*(y4 - y6)^2)/denom_elas2;
    a(4,5)=(lambda*(x3 - x4)*(y4 - y6))/denom_elas2;
    a(4,6)=(-2*(lambda + mu)*(x3 - x4)*(x4 - x6) - (lambda + 2*mu)*(y3 - y4)*(y4 - y6))/denom_elas2;
%     a(5,1)=0;
%     a(5,2)=0;
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/denom_elas2;
    a(5,6)=-((lambda*(x3 - x4)*(y3 - y4))/denom_elas2);
%     a(6,1)=0;
%     a(6,2)=0;
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/denom_elas2;
    
    b=(1/2)*[-1 0 0 0;0 -1 0 0;1 0 -1 0;0 1 0 -1;0 0 -1 0;0 0 0 -1];
    
    c=(1/4)*[y3-y4 x4-x3 y6-y1 x1-x6 y3-y4 x4-x3]';
    
    d=zeros(vdim,2);
    d(1,1)=(lambda + 2*mu)*(x3 - x4)^2 + lambda*(x3 - x4)*(y1 - y4) + 2*(lambda + mu)*(y3 - y4)^2;
    d(2,1)=-2*(lambda + mu)*(x1 - x4)*(x3 - x4) - lambda*(x3 - x4)*(y3 - y4) - (lambda + 2*mu)*(y1 - y4)*(y3 - y4);
    d(3,1)=(lambda + 2*mu)*(x1 - x4)*(-x3 + x4) + lambda*(x1 - x4)*(-y1 + y4) + 2*(lambda + mu)*(y1 - y4)*(-y3 + y4); 
    d(4,1)=2*(lambda + mu)*(x1 - x4)^2 + lambda*(x3 - x4)*(y1 - y4) + (lambda + 2*mu)*(y1 - y4)^2;
%     d(5,1)=0;
%     d(6,1)=0;
%     d(1,2)=0;
%     d(2,2)=0;
    d(3,2)=-((lambda + 2*mu)*(x3 - x4)*(x4 - x6)) + lambda*(-x4 + x6)*(y4 - y6) - 2*(lambda + mu)*(y3 - y4)*(y4 - y6); 
    d(4,2)=2*(lambda + mu)*(x4 - x6)^2 + lambda*(x3 - x4)*(y4 - y6) + (lambda + 2*mu)*(y4 - y6)^2;
    d(5,2)=(lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2 + lambda*(x3 - x4)*(y4 - y6);
    d(6,2)=-2*(lambda + mu)*(x3 - x4)*(x4 - x6) - lambda*(x3 - x4)*(y3 - y4) - (lambda + 2*mu)*(y3 - y4)*(y4 - y6);
    
    d(:,1)=denom_d1*d(:,1);
    d(:,2)=denom_d2*d(:,2);
    
    % CC.D en el nodo Oeste, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
        % frontera Oeste 1
            Pg1u_L4=int_dir_simpson(i,j,i,j-1,1);
            Pg2u_L4=int_dir_simpson(i,j,i,j-1,2);
        % frontera Oeste 2
            Pg1u_L11=int_dir_simpson(i,j,i,j+1,1);
            Pg2u_L11=int_dir_simpson(i,j,i,j+1,2);
        Pgu=(1/2)*[-Pg1u_L4;-Pg2u_L4;0;0;-Pg1u_L11;-Pg2u_L11];
   
    pp(1,1)=p(i,j-1);
    pp(2,1)=p(i,j);
    gg=gamma(i,j);
    
    ss=a\(Pgu-b*u([ind1u:ind2u,ind3u:ind4u],1)-c*gg-d*pp);
    sigma(ld+1:ld+vdim)=ss;
    ld=ld+vdim;
    
    sx1(i,j,1)=ss(1,1);
    sx2(i,j,1)=ss(2,1);
    sy1(i,j,1)=ss(3,1);
    sy2(i,j,1)=ss(4,1);
    sx1(i,j,2)=ss(5,1);
    sx2(i,j,2)=ss(6,1);
    
    % Central nodes 
    vdim=8;
    pdim=4;
    pp=zeros(pdim,1);

    a=zeros(vdim,vdim);
    d=zeros(vdim,4);
    % Matrix b is the same for every central node
    b=(1/2)*[1 0 -1 0 0 0 0 0;0 1 0 -1 0 0 0 0;0 0 1 0 0 0 -1 0;...
             0 0 0 1 0 0 0 -1;0 0 0 0 1 0 -1 0;0 0 0 0 0 1 0 -1;...
             1 0 0 0 -1 0 0 0;0 1 0 0 0 -1 0 0;];
         
    for i=2:N
        ind2u=(i-1+(j-2)*N)*2;
        ind1u=ind2u-1;
        ind4u=(i+(j-2)*N)*2;
%         ind3u=ind4u-1;
        ind6u=(i-1+(j-1)*N)*2;
        ind5u=ind6u-1;
        ind8u=(i+(j-1)*N)*2;
%         ind7u=ind8u-1;
        
%         x1=x(i-1,j-1);
%         y1=y(i-1,j-1);
        x2=x(i,j-1);
    y2=y(i,j-1);
    x3=x(i,j);
    y3=y(i,j);
    x4=x(i-1,j);
    y4=y(i-1,j);
%     x5=x(i+1,j-1);
%     y5=y(i+1,j-1);
    x6=x(i+1,j);
    y6=y(i+1,j);
%     x7=x(i+1,j+1);
%     y7=y(i+1,j+1);
    x8=x(i,j+1);
    y8=y(i,j+1);
%     x9=x(i-1,j+1);
%     y9=y(i-1,j+1);
    
    JE1r3=abs(x4*(y2 - y3) + x2*(y3 - y4) + x3*(-y2 + y4));
    denom_elas1=coef_denom*JE1r3;
    denom_d1=coef_d/JE1r3;
    JE2r4=abs(x6*(-y2 + y3) + x3*(y2 - y6) + x2*(-y3 + y6));
    denom_elas2=coef_denom*JE2r4;
    denom_d2=coef_d/JE2r4;
    JE3r1=abs(x8*(y3 - y6) + x3*(y6 - y8) + x6*(-y3 + y8));
    denom_elas3=coef_denom*JE3r1;
    denom_d3=coef_d/JE3r1;
    JE4r2=abs(x8*(-y3 + y4) + x4*(y3 - y8) + x3*(-y4 + y8));
    denom_elas4=coef_denom*JE4r2;
    denom_d4=coef_d/JE4r2;
    
    a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/denom_elas1 + ...
           ((lambda + 2*mu)*(x3 - x6)^2 + 2*(lambda + mu)*(y3 - y6)^2)/denom_elas2;
    a(1,2)=-((lambda*(x3 - x4)*(y3 - y4))/denom_elas1) - (lambda*(x3 - x6)*(y3 - y6))/denom_elas2;
    a(1,3)=((lambda + 2*mu)*(x2 - x3)*(x3 - x6) + 2*(lambda + mu)*(y2 - y3)*(y3 - y6))/denom_elas2;
    a(1,4)=-((lambda*(x3 - x6)*(y2 - y3))/denom_elas2);
%     a(1,5)=0;
%     a(1,6)=0;
    a(1,7)=-(((lambda + 2*mu)*(x2 - x3)*(x3 - x4) + 2*(lambda + mu)*(y2 - y3)*(y3 - y4))/denom_elas1);
    a(1,8)=(lambda*(x3 - x4)*(y2 - y3))/denom_elas1;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/denom_elas1 + ...
           (2*(lambda + mu)*(x3 - x6)^2 + (lambda + 2*mu)*(y3 - y6)^2)/denom_elas2;
    a(2,3)=-((lambda*(x2 - x3)*(y3 - y6))/denom_elas2);
    a(2,4)=(2*(lambda + mu)*(x2 - x3)*(x3 - x6) + (lambda + 2*mu)*(y2 - y3)*(y3 - y6))/denom_elas2;
%     a(2,5)=0;
%     a(2,6)=0;
    a(2,7)=(lambda*(x2 - x3)*(y3 - y4))/denom_elas1;
    a(2,8)=(-2*(lambda + mu)*(x2 - x3)*(x3 - x4) - (lambda + 2*mu)*(y2 - y3)*(y3 - y4))/denom_elas1;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/denom_elas2 + ...
           ((lambda + 2*mu)*(x3 - x8)^2 + 2*(lambda + mu)*(y3 - y8)^2)/denom_elas3;
    a(3,4)=-((lambda*(x2 - x3)*(y2 - y3))/denom_elas2) - (lambda*(x3 - x8)*(y3 - y8))/denom_elas3;
    a(3,5)=((lambda + 2*mu)*(x3 - x6)*(x3 - x8) + 2*(lambda + mu)*(y3 - y6)*(y3 - y8))/denom_elas3;
    a(3,6)=-((lambda*(x3 - x8)*(y3 - y6))/denom_elas3);
%     a(3,7)=0;
%     a(3,8)=0;
    a(4,1)=a(1,4);
    a(4,2)=a(2,4);
    a(4,3)=a(3,4);
    a(4,4)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/denom_elas2 + ...
           (2*(lambda + mu)*(x3 - x8)^2 + (lambda + 2*mu)*(y3 - y8)^2)/denom_elas3;
    a(4,5)=-((lambda*(x3 - x6)*(y3 - y8))/denom_elas3);
    a(4,6)=(2*(lambda + mu)*(x3 - x6)*(x3 - x8) + (lambda + 2*mu)*(y3 - y6)*(y3 - y8))/denom_elas3;
%     a(4,7)=0;
%     a(4,8)=0;
%     a(5,1)=0;
%     a(5,2)=0;
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/denom_elas4 + ...
           ((lambda + 2*mu)*(x3 - x6)^2 + 2*(lambda + mu)*(y3 - y6)^2)/denom_elas3;
    a(5,6)=-((lambda*(x3 - x4)*(y3 - y4))/denom_elas4) - (lambda*(x3 - x6)*(y3 - y6))/denom_elas3;
    a(5,7)=-(((lambda + 2*mu)*(x3 - x4)*(x3 - x8) + 2*(lambda + mu)*(y3 - y4)*(y3 - y8))/denom_elas4);
    a(5,8)=(lambda*(x3 - x4)*(y3 - y8))/denom_elas4;
%     a(6,1)=0;
%     a(6,2)=0;
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/denom_elas4 + ...
           (2*(lambda + mu)*(x3 - x6)^2 + (lambda + 2*mu)*(y3 - y6)^2)/denom_elas3;
    a(6,7)=(lambda*(x3 - x8)*(y3 - y4))/denom_elas4;
    a(6,8)=-((2*(lambda + mu)*(x3 - x4)*(x3 - x8) + (lambda + 2*mu)*(y3 - y4)*(y3 - y8))/denom_elas4);
    a(7,1)=a(1,7);
    a(7,2)=a(2,7);
%     a(7,3)=0;
%     a(7,4)=0;
    a(7,5)=a(5,7);
    a(7,6)=a(6,7);
    a(7,7)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/denom_elas1 + ...
           ((lambda + 2*mu)*(x3 - x8)^2 + 2*(lambda + mu)*(y3 - y8)^2)/denom_elas4;
    a(7,8)=-((lambda*(x2 - x3)*(y2 - y3))/denom_elas1) - (lambda*(x3 - x8)*(y3 - y8))/denom_elas4;
    a(8,1)=a(1,8);
    a(8,2)=a(2,8);
%     a(8,3)=0;
%     a(8,4)=0;
    a(8,5)=a(5,8);
    a(8,6)=a(6,8);
    a(8,7)=a(7,8);
    a(8,8)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/denom_elas1 + ...
           (2*(lambda + mu)*(x3 - x8)^2 + (lambda + 2*mu)*(y3 - y8)^2)/denom_elas4;
    
    c=(1/4)*[y6-y4 x4-x6 y8-y2 x2-x8 y6-y4 x4-x6 y8-y2 x2-x8]';
    
    % ¡¡OJO con las posiciones 3 y 4 (van intercaladas <- va antes la 4 que 
    % la 3!!
    d(1,1)=(lambda + 2*mu)*(x3 - x4)^2 + lambda*(x3 - x4)*(y2 - y3) + 2*(lambda + mu)*(y3 - y4)^2;
    d(2,1)=-2*(lambda + mu)*(x2 - x3)*(x3 - x4) + lambda*(-x3 + x4)*(y3 - y4) + (lambda + 2*mu)*(-y2 + y3)*(y3 - y4);
%     d(3,1)=0;
%     d(4,1)=0;
%     d(5,1)=0;
%     d(6,1)=0;
    d(7,1)=-((lambda + 2*mu)*(x2 - x3)*(x3 - x4)) - lambda*(x2 - x3)*(y2 - y3) - 2*(lambda + mu)*(y2 - y3)*(y3 - y4);
    d(8,1)=2*(lambda + mu)*(x2 - x3)^2 + lambda*(x3 - x4)*(y2 - y3) + (lambda + 2*mu)*(y2 - y3)^2;
    d(1,2)=(lambda + 2*mu)*(x3 - x6)^2 + lambda*(-x3 + x6)*(y2 - y3) + 2*(lambda + mu)*(y3 - y6)^2;
    d(2,2)=2*(lambda + mu)*(x2 - x3)*(x3 - x6) + lambda*(-x3 + x6)*(y3 - y6) + (lambda + 2*mu)*(y2 - y3)*(y3 - y6);
    d(3,2)=(lambda + 2*mu)*(x2 - x3)*(x3 - x6) + lambda*(-x2 + x3)*(y2 - y3) + 2*(lambda + mu)*(y2 - y3)*(y3 - y6);
    d(4,2)=2*(lambda + mu)*(x2 - x3)^2 + lambda*(-x3 + x6)*(y2 - y3) + (lambda + 2*mu)*(y2 - y3)^2;
%     d(5,2)=0;
%     d(6,2)=0;
%     d(7,2)=0;
%     d(8,2)=0;

%     ¡OJO!: ¡Posicones 3 y 4 van intercaladas! (P4 va antes que P3)

%     d(1,3)=0;
%     d(2,3)=0;
%     d(3,3)=0;
%     d(4,3)=0;
    d(5,3)=(lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2 + lambda*(x3 - x4)*(y3 - y8);
    d(6,3)=-2*(lambda + mu)*(x3 - x4)*(x3 - x8) - lambda*(x3 - x4)*(y3 - y4) - (lambda + 2*mu)*(y3 - y4)*(y3 - y8);
    d(7,3)=-((lambda + 2*mu)*(x3 - x4)*(x3 - x8)) - lambda*(x3 - x8)*(y3 - y8) - 2*(lambda + mu)*(y3 - y4)*(y3 - y8);
    d(8,3)=2*(lambda + mu)*(x3 - x8)^2 + lambda*(x3 - x4)*(y3 - y8) + (lambda + 2*mu)*(y3 - y8)^2;

%     d(1,4)=0;
%     d(2,4)=0;
    d(3,4)=(lambda + 2*mu)*(x3 - x6)*(x3 - x8) + lambda*(-x3 + x8)*(y3 - y8) + 2*(lambda + mu)*(y3 - y6)*(y3 - y8);
    d(4,4)=2*(lambda + mu)*(x3 - x8)^2 + lambda*(-x3 + x6)*(y3 - y8) + (lambda + 2*mu)*(y3 - y8)^2;
    d(5,4)=(lambda + 2*mu)*(x3 - x6)^2 + 2*(lambda + mu)*(y3 - y6)^2 + lambda*(-x3 + x6)*(y3 - y8);
    d(6,4)=2*(lambda + mu)*(x3 - x6)*(x3 - x8) + lambda*(-x3 + x6)*(y3 - y6) + (lambda + 2*mu)*(y3 - y6)*(y3 - y8);
%     d(7,4)=0;
%     d(8,4)=0;
    
    d(:,1)=denom_d1*d(:,1);
    d(:,2)=denom_d2*d(:,2);    
%     ¡OJO!: Índices 3 y 4 van intercalados    
    d(:,3)=denom_d4*d(:,3);
    d(:,4)=denom_d3*d(:,4);

        % NO Dir.BC
        pp(1,1)=p(i-1,j-1);
        pp(2,1)=p(i,j-1);
        pp(3,1)=p(i-1,j);
        pp(4,1)=p(i,j);
        gg=gamma(i,j);

        ss=a\(-b*u([ind1u:ind4u,ind5u:ind8u],1)-c*gg-d*pp);
        sigma(ld+1:ld+vdim)=ss;
        ld=ld+vdim;

        sx1(i,j,1)=ss(1,1);
        sx2(i,j,1)=ss(2,1);
        sy1(i,j,1)=ss(3,1);
        sy2(i,j,1)=ss(4,1);
        sx1(i,j,2)=ss(5,1);
        sx2(i,j,2)=ss(6,1);
        sy1(i,j,2)=ss(7,1);
        sy2(i,j,2)=ss(8,1);
    end
    
    % East nodes
    i=N+1;
    ind2u=(i-1+(j-2)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i-1+(j-1)*N)*2;
    ind3u=ind4u-1;
    
    vdim=6;
    pdim=2;
    pp=zeros(pdim,1);
    
%     x1=x(i-1,j-1);
%     y1=y(i-1,j-1);
    x2=x(i,j-1);
    y2=y(i,j-1);
    x3=x(i,j);
    y3=y(i,j);
    x4=x(i-1,j);
    y4=y(i-1,j);
    x5=x(i,j+1);
    y5=y(i,j+1);
%     x6=x(i-1,j+1);
%     y6=y(i-1,j+1);   
    
    JE1r3=abs(x4*(y2 - y3) + x2*(y3 - y4) + x3*(-y2 + y4));
    denom_elas1=coef_denom*JE1r3;
    denom_d1=coef_d/JE1r3;
    JE2r2=abs(x5*(-y3 + y4) + x4*(y3 - y5) + x3*(-y4 + y5));
    denom_elas2=coef_denom*JE2r2;
    denom_d2=coef_d/JE2r2;
    
    a=zeros(vdim,vdim);
    a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/denom_elas1;
    a(1,2)=-((lambda*(x3 - x4)*(y3 - y4))/denom_elas1);
%     a(1,3)=0;
%     a(1,4)=0;
    a(1,5)=-(((lambda + 2*mu)*(x2 - x3)*(x3 - x4) + 2*(lambda + mu)*(y2 - y3)*(y3 - y4))/denom_elas1);
    a(1,6)=(lambda*(x3 - x4)*(y2 - y3))/denom_elas1;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/denom_elas1;
%     a(2,3)=0;
%     a(2,4)=0:
    a(2,5)=(lambda*(x2 - x3)*(y3 - y4))/denom_elas1;
    a(2,6)=-((2*(lambda + mu)*(x2 - x3)*(x3 - x4) + (lambda + 2*mu)*(y2 - y3)*(y3 - y4))/denom_elas1);
%     a(3,1)=0;
%     a(3,2)=0;
    a(3,3)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/denom_elas2;
    a(3,4)=-((lambda*(x3 - x4)*(y3 - y4))/denom_elas2);
    a(3,5)=-(((lambda + 2*mu)*(x3 - x4)*(x3 - x5) + 2*(lambda + mu)*(y3 - y4)*(y3 - y5))/denom_elas2);
    a(3,6)=(lambda*(x3 - x4)*(y3 - y5))/denom_elas2;
%     a(4,1)=0;
%     a(4,2)=0;
    a(4,3)=a(3,4);
    a(4,4)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/denom_elas2;
    a(4,5)=(lambda*(x3 - x5)*(y3 - y4))/denom_elas2;
    a(4,6)=-((2*(lambda + mu)*(x3 - x4)*(x3 - x5) + (lambda + 2*mu)*(y3 - y4)*(y3 - y5))/denom_elas2);
    a(5,1)=a(1,5);
    a(5,2)=a(2,5);
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/denom_elas1 + ...
           ((lambda + 2*mu)*(x3 - x5)^2 + 2*(lambda + mu)*(y3 - y5)^2)/denom_elas2;
    a(5,6)=-((lambda*(x2 - x3)*(y2 - y3))/denom_elas1) - (lambda*(x3 - x5)*(y3 - y5))/denom_elas2;
    a(6,1)=a(1,6);
    a(6,2)=a(2,6);
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/denom_elas1 +...
           (2*(lambda + mu)*(x3 - x5)^2 + (lambda + 2*mu)*(y3 - y5)^2)/denom_elas2;
    
    b=(1/2)*[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1;1 0 -1 0;0 1 0 -1];
    
    c=(1/4)*[y3-y4 x4-x3 y3-y4 x4-x3 y5-y2 x2-x5]';
    
    d=zeros(vdim,2);
    d(1,1)=(lambda + 2*mu)*(x3 - x4)^2 + lambda*(x3 - x4)*(y2 - y3) + 2*(lambda + mu)*(y3 - y4)^2;
    d(2,1)=-2*(lambda + mu)*(x2 - x3)*(x3 - x4) + lambda*(-x3 + x4)*(y3 - y4) + (lambda + 2*mu)*(-y2 + y3)*(y3 - y4);
%     d(3,1)=0;
%     d(4,1)=0;
    d(5,1)=-((lambda + 2*mu)*(x2 - x3)*(x3 - x4)) - lambda*(x2 - x3)*(y2 - y3) - 2*(lambda + mu)*(y2 - y3)*(y3 - y4);
    d(6,1)=2*(lambda + mu)*(x2 - x3)^2 + lambda*(x3 - x4)*(y2 - y3) + (lambda + 2*mu)*(y2 - y3)^2;
%     d(1,2)=0;
%     d(2,2)=0;
    d(3,2)=(lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2 + lambda*(x3 - x4)*(y3 - y5);
    d(4,2)=-2*(lambda + mu)*(x3 - x4)*(x3 - x5) - lambda*(x3 - x4)*(y3 - y4) - (lambda + 2*mu)*(y3 - y4)*(y3 - y5);
    d(5,2)=-((lambda + 2*mu)*(x3 - x4)*(x3 - x5)) + lambda*(-x3 + x5)*(y3 - y5) - 2*(lambda + mu)*(y3 - y4)*(y3 - y5);
    d(6,2)=2*(lambda + mu)*(x3 - x5)^2 + lambda*(x3 - x4)*(y3 - y5) + (lambda + 2*mu)*(y3 - y5)^2;
    
    d(:,1)=denom_d1*d(:,1);
    d(:,2)=denom_d2*d(:,2);
    
    % CC.D en el nodo Este, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
        % frontera Este 1
            Pg1u_L7=int_dir_simpson(i,j,i,j-1,1);
            Pg2u_L7=int_dir_simpson(i,j,i,j-1,2);
        % frontera Este 2
            Pg1u_L14=int_dir_simpson(i,j,i,j+1,1);
            Pg2u_L14=int_dir_simpson(i,j,i,j+1,2);
        Pgu=(1/2)*[Pg1u_L7;Pg2u_L7;Pg1u_L14;Pg2u_L14;0;0];

    pp(1,1)=p(i-1,j-1);
    pp(2,1)=p(i-1,j);
    gg=gamma(i,j);
    
    ss=a\(Pgu-b*u([ind1u:ind2u,ind3u:ind4u],1)-c*gg-d*pp);
    sigma(ld+1:ld+vdim)=ss;
    ld=ld+vdim;
    
    sx1(i,j,1)=ss(1,1);
    sx2(i,j,1)=ss(2,1);
    sx1(i,j,2)=ss(3,1);
    sx2(i,j,2)=ss(4,1);
    sy1(i,j,2)=ss(5,1);
    sy2(i,j,2)=ss(6,1);
end

% North-West corner node
i=1;
j=N+1;
ind2u=(i+(j-2)*N)*2;
ind1u=ind2u-1;
vdim=4;
pdim=1;
pp=zeros(pdim,1);

x1=x(i,j-1);
y1=y(i,j-1);           
% x2=x(i+1,j-1);
% y2=y(i+1,j-1);
x3=x(i+1,j);
y3=y(i+1,j);
x4=x(i,j);
y4=y(i,j);

Jer4=abs(x4*(y1 - y3) + x1*(y3 - y4) + x3*(-y1 + y4));
denom_elas=coef_denom*Jer4;
denom_d=coef_d/Jer4;

a=zeros(vdim,vdim);
a(1,1)=(lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2;
a(1,2)=lambda*(-x3 + x4)*(y3 - y4);
a(1,3)=-((lambda + 2*mu)*(x1 - x4)*(x3 - x4)) - 2*(lambda + mu)*(y1 - y4)*(y3 - y4);
a(1,4)=lambda*(x3 - x4)*(y1 - y4);
a(2,1)=a(1,2);
a(2,2)=2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2;
a(2,3)=lambda*(x1 - x4)*(y3 - y4);
a(2,4)=-2*(lambda + mu)*(x1 - x4)*(x3 - x4) - (lambda + 2*mu)*(y1 - y4)*(y3 - y4);
a(3,1)=a(1,3);
a(3,2)=a(2,3);
a(3,3)=(lambda + 2*mu)*(x1 - x4)^2 + 2*(lambda + mu)*(y1 - y4)^2;
a(3,4)=lambda*(-x1 + x4)*(y1 - y4);
a(4,1)=a(1,4);
a(4,2)=a(2,4);
a(4,3)=a(3,4);
a(4,4)=2*(lambda + mu)*(x1 - x4)^2 + (lambda + 2*mu)*(y1 - y4)^2;
a=(1/denom_elas)*a;

b=(1/2)*[-1 0;0 -1;1 0;0 1];

c=(1/4)*[y3-y4 x4-x3 y4-y1 x1-x4]';

d=zeros(vdim,1);
d(1,1)=(lambda + 2*mu)*(x3 - x4)^2 + lambda*(x3 - x4)*(y1 - y4) + 2*(lambda + mu)*(y3 - y4)^2;
d(2,1)=-2*(lambda + mu)*(x1 - x4)*(x3 - x4) - lambda*(x3 - x4)*(y3 - y4) - (lambda + 2*mu)*(y1 - y4)*(y3 - y4);
d(3,1)=-((lambda + 2*mu)*(x1 - x4)*(x3 - x4)) - lambda*(x1 - x4)*(y1 - y4) - 2*(lambda + mu)*(y1 - y4)*(y3 - y4);
d(4,1)=2*(lambda + mu)*(x1 - x4)^2 + lambda*(x3 - x4)*(y1 - y4) + (lambda + 2*mu)*(y1 - y4)^2;
d=denom_d*d;

% CC.D en el nodo N-W, fórmula de cuadratura para los término Pg1 y Pg2 en:
    % frontera Oeste
        Pg1u_L18=int_dir_simpson(i,j,i,j-1,1);
        Pg2u_L18=int_dir_simpson(i,j,i,j-1,2);
    % frontera Norte
        Pg1u_L22=int_dir_simpson(i,j,i+1,j,1);
        Pg2u_L22=int_dir_simpson(i,j,i+1,j,2);
    Pgu=(1/2)*[-Pg1u_L18;-Pg2u_L18;Pg1u_L22;Pg2u_L22];
    
pp(1,1)=p(i,j-1);
gg=gamma(i,j);

ss=a\(Pgu-b*u(ind1u:ind2u)-c*gg-d*pp);
sigma(ld+1:ld+vdim)=ss;
ld=ld+vdim;

sx1(i,j,1)=ss(1,1);
sx2(i,j,1)=ss(2,1);
sy1(i,j,1)=ss(3,1);
sy2(i,j,1)=ss(4,1);

% North nodes (j=N+1)
vdim=6;
pdim=2;
pp=zeros(pdim,1);

% Matrix b is the same for every North node
b=(1/2)*[1 0 -1 0;0 1 0 -1;0 0 1 0;0 0 0 1;1 0 0 0;0 1 0 0];

a=zeros(vdim,vdim);
d=zeros(vdim,2);

for i=2:N
    ind2u=(i-1+(j-2)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i+(j-2)*N)*2;
%     ind3u=ind4u-1;
    
%     x1=x(i-1,j-1);
%     y1=y(i-1,j-1);
    x2=x(i,j-1);
    y2=y(i,j-1);
    x3=x(i,j);
    y3=y(i,j);
    x4=x(i-1,j);
    y4=y(i-1,j);
%     x5=x(i+1,j-1);
%     y5=y(i+1,j-1);
    x6=x(i+1,j);
    y6=y(i+1,j);   
    
    JE1r3=abs(x4*(y2 - y3) + x2*(y3 - y4) + x3*(-y2 + y4));
    denom_elas1=coef_denom*JE1r3;
    denom_d1=coef_d/JE1r3;
    JE2r4=abs(x6*(-y2 + y3) + x3*(y2 - y6) + x2*(-y3 + y6));
    denom_elas2=coef_denom*JE2r4;
    denom_d2=coef_d/JE2r4;
    
    a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/denom_elas1 + ...
           ((lambda + 2*mu)*(x3 - x6)^2 + 2*(lambda + mu)*(y3 - y6)^2)/denom_elas2;
    a(1,2)=-((lambda*(x3 - x4)*(y3 - y4))/denom_elas1) - (lambda*(x3 - x6)*(y3 - y6))/denom_elas2;
    a(1,3)=((lambda + 2*mu)*(x2 - x3)*(x3 - x6) + 2*(lambda + mu)*(y2 - y3)*(y3 - y6))/denom_elas2;
    a(1,4)=-((lambda*(x3 - x6)*(y2 - y3))/denom_elas2);
    a(1,5)=-(((lambda + 2*mu)*(x2 - x3)*(x3 - x4) + 2*(lambda + mu)*(y2 - y3)*(y3 - y4))/denom_elas1);
    a(1,6)=(lambda*(x3 - x4)*(y2 - y3))/denom_elas1;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/denom_elas1 + ...
           (2*(lambda + mu)*(x3 - x6)^2 + (lambda + 2*mu)*(y3 - y6)^2)/denom_elas2;
    a(2,3)=-((lambda*(x2 - x3)*(y3 - y6))/denom_elas2);
    a(2,4)=(2*(lambda + mu)*(x2 - x3)*(x3 - x6) + (lambda + 2*mu)*(y2 - y3)*(y3 - y6))/denom_elas2;
    a(2,5)=(lambda*(x2 - x3)*(y3 - y4))/denom_elas1;
    a(2,6)=(-2*(lambda + mu)*(x2 - x3)*(x3 - x4) - (lambda + 2*mu)*(y2 - y3)*(y3 - y4))/denom_elas1;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/denom_elas2;
    a(3,4)=-((lambda*(x2 - x3)*(y2 - y3))/denom_elas2);
%     a(3,5)=0;
%     a(3,6)=0;
    a(4,1)= a(1,4);
    a(4,2)= a(2,4);
    a(4,3)= a(3,4);
    a(4,4)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/denom_elas2;
%     a(4,5)=0;
%     a(4,6)=0;
    a(5,1)=a(1,5);
    a(5,2)=a(2,5);
%     a(5,3)=0;
%     a(5,4)=0;
    a(5,5)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/denom_elas1;
    a(5,6)=-((lambda*(x2 - x3)*(y2 - y3))/denom_elas1);
    a(6,1)=a(1,6);
    a(6,2)=a(2,6);
%     a(6,3)=0;
%     a(6,4)=0;
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/denom_elas1;
    
    c=(1/4)*[y6-y4 x4-x6 y3-y2 x2-x3 y3-y2 x2-x3]';
    
    d(1,1)=(lambda + 2*mu)*(x3 - x4)^2 + lambda*(x3 - x4)*(y2 - y3) + 2*(lambda + mu)*(y3 - y4)^2;
    d(2,1)=-2*(lambda + mu)*(x2 - x3)*(x3 - x4) + lambda*(-x3 + x4)*(y3 - y4) + (lambda + 2*mu)*(-y2 + y3)*(y3 - y4);
%     d(3,1)=0;
%     d(4,1)=0;
    d(5,1)=-((lambda + 2*mu)*(x2 - x3)*(x3 - x4)) - lambda*(x2 - x3)*(y2 - y3) - 2*(lambda + mu)*(y2 - y3)*(y3 - y4);
    d(6,1)=2*(lambda + mu)*(x2 - x3)^2 + lambda*(x3 - x4)*(y2 - y3) + (lambda + 2*mu)*(y2 - y3)^2;
    d(1,2)=(lambda + 2*mu)*(x3 - x6)^2 + lambda*(-x3 + x6)*(y2 - y3) + 2*(lambda + mu)*(y3 - y6)^2;
    d(2,2)=2*(lambda + mu)*(x2 - x3)*(x3 - x6) + lambda*(-x3 + x6)*(y3 - y6) + (lambda + 2*mu)*(y2 - y3)*(y3 - y6);
    d(3,2)=(lambda + 2*mu)*(x2 - x3)*(x3 - x6) + lambda*(-x2 + x3)*(y2 - y3) + 2*(lambda + mu)*(y2 - y3)*(y3 - y6);
    d(4,2)=2*(lambda + mu)*(x2 - x3)^2 + lambda*(-x3 + x6)*(y2 - y3) + (lambda + 2*mu)*(y2 - y3)^2;
%     d(5,2)=0;
%     d(6,2)=0;

    d(:,1)=denom_d1*d(:,1);
    d(:,2)=denom_d2*d(:,2);
    
    % CC.D en el nodo Norte, fórmula de cuadratura para los términos Pg1 y Pg2 de u en:
        % frontera Norte 1
            Pg1u_L23=int_dir_simpson(i,j,i+1,j,1);
            Pg2u_L23=int_dir_simpson(i,j,i+1,j,2);
        % frontera Norte 2
            Pg1u_L22=int_dir_simpson(i-1,j,i,j,1);
            Pg2u_L22=int_dir_simpson(i-1,j,i,j,2);
        Pgu=(1/2)*[0;0;Pg1u_L23;Pg2u_L23;Pg1u_L22;Pg2u_L22];
        
        pp(1,1)=p(i-1,j-1);
        pp(2,1)=p(i,j-1);
        gg=gamma(i,j);
        
        ss=a\(Pgu-b*u(ind1u:ind4u)-c*gg-d*pp);
        sigma(ld+1:ld+vdim)=ss;
        ld=ld+vdim;
        
        sx1(i,j,1)=ss(1,1);
        sx2(i,j,1)=ss(2,1);
        sy1(i,j,1)=ss(3,1);
        sy2(i,j,1)=ss(4,1);
        sy1(i,j,2)=ss(5,1);
        sy2(i,j,2)=ss(6,1);
end

% North-East corner node (j=N+1)
i=N+1;
ind2u=(i-1+(j-2)*N)*2;
ind1u=ind2u-1;
vdim=4;
pdim=1;
pp=zeros(pdim,1);

% x1=x(i-1,j-1);  
% y1=y(i-1,j-1);            
x2=x(i,j-1);
y2=y(i,j-1);
x3=x(i,j);
y3=y(i,j);
x4=x(i-1,j);
y4=y(i-1,j);

Jer3=abs(x4*(y2 - y3) + x2*(y3 - y4) + x3*(-y2 + y4));
denom_elas=coef_denom*Jer3;
denom_d=coef_d/Jer3;

a=zeros(vdim,vdim);
a(1,1)=(lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2;
a(1,2)=-lambda*(x3 - x4)*(y3 - y4);
a(1,3)=-((lambda + 2*mu)*(x2 - x3)*(x3 - x4) + 2*(lambda + mu)*(y2 - y3)*(y3 - y4));
a(1,4)=lambda*(x3 - x4)*(y2 - y3);
a(2,1)=a(1,2);
a(2,2)=2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2;
a(2,3)=lambda*(x2 - x3)*(y3 - y4);
a(2,4)=-(2*(lambda + mu)*(x2 - x3)*(x3 - x4) + (lambda + 2*mu)*(y2 - y3)*(y3 - y4));
a(3,1)=a(1,3);
a(3,2)=a(2,3);
a(3,3)=(lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2;
a(3,4)=-lambda*(x2 - x3)*(y2 - y3);
a(4,1)=a(1,4);
a(4,2)=a(2,4);
a(4,3)=a(3,4);
a(4,4)=2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2;
a=(1/denom_elas)*a;

b=(1/2)*[1 0;0 1;1 0;0 1];

c=(1/4)*[y3-y4 x4-x3 y3-y2 x2-x3]';

d=zeros(vdim,1);
d(1,1)=(lambda + 2*mu)*(x3 - x4)^2 + lambda*(x3 - x4)*(y2 - y3) + 2*(lambda + mu)*(y3 - y4)^2;
d(2,1)=-2*(lambda + mu)*(x2 - x3)*(x3 - x4) + lambda*(-x3 + x4)*(y3 - y4) + (lambda + 2*mu)*(-y2 + y3)*(y3 - y4);
d(3,1)=-((lambda + 2*mu)*(x2 - x3)*(x3 - x4)) - lambda*(x2 - x3)*(y2 - y3) - 2*(lambda + mu)*(y2 - y3)*(y3 - y4);
d(4,1)=2*(lambda + mu)*(x2 - x3)^2 + lambda*(x3 - x4)*(y2 - y3) + (lambda + 2*mu)*(y2 - y3)^2;
d=denom_d*d;

% CC.D en el nodo N-E, fórmula de cuadratura para los términos Pg1 y Pg2 de u en:
    % frontera Este
        Pg1u_L21=int_dir_simpson(i,j,i,j-1,1);
        Pg2u_L21=int_dir_simpson(i,j,i,j-1,2);
    % frontera Norte
        Pg1u_L24=int_dir_simpson(i-1,j,i,j,1);
        Pg2u_L24=int_dir_simpson(i-1,j,i,j,2);
    Pgu=(1/2)*[Pg1u_L21;Pg2u_L21;Pg1u_L24;Pg2u_L24];

pp(1,1)=p(i-1,j-1);
gg=gamma(i,j);

ss=a\(Pgu-b*u(ind1u:ind2u)-c*gg-d*pp);
sigma(ld+1:ld+vdim)=ss;

sx1(i,j,1)=ss(1,1);
sx2(i,j,1)=ss(2,1);
sy1(i,j,2)=ss(3,1);
sy2(i,j,2)=ss(4,1);

return
end
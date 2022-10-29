function gamma=compute_gamma(u,p,tt)

global NN x y alpha nu 

N=NN;

gamma=zeros((N+1)*(N+1),1);

coef_d=(alpha*(-1 + nu + 2*nu^2))/20;

% South-West corner node
i=1;
j=1;
ind2u=(i+(j-1)*N)*2;
ind1u=ind2u-1;
ind1p=ind2u/2;
vdim=4;
pos=1;

x1=x(i,j);  
y1=y(i,j);            
x2=x(i+1,j);
y2=y(i+1,j);
% x3=x(i+1,j+1);
% y3=y(i+1,j+1);
x4=x(i,j+1);
y4=y(i,j+1);

Jer1=abs(x4*(y1 - y2) + x1*(y2 - y4) + x2*(-y1 + y4));
denom=20*Jer1;
coef_a=-1/denom;

% Local matrix (A sigma, sigma)
a=zeros(vdim,vdim);
a(1,1)=(1 + nu)*((-1 + nu)*(x1 - x4)^2 - (y1 - y4)^2);
a(1,2)=nu*(1 + nu)*(x1 - x4)*(y1 - y4);
a(1,3)=(1 + nu)*((-1 + nu)*(x1 - x2)*(x1 - x4) - (y1 - y2)*(y1 - y4));
a(1,4)=nu*(1 + nu)*(x1 - x4)*(y1 - y2);
a(2,1)=a(1,2);
a(2,2)=(1 + nu)*(-(x1 - x4)^2 + (-1 + nu)*(y1 - y4)^2);
a(2,3)=nu*(1 + nu)*(x1 - x2)*(y1 - y4);
a(2,4)=(1 + nu)*(-((x1 - x2)*(x1 - x4)) + (-1 + nu)*(y1 - y2)*(y1 - y4));
a(3,1)=a(1,3);
a(3,2)=a(2,3);
a(3,3)=(1 + nu)*((-1 + nu)*(x1 - x2)^2 - (y1 - y2)^2);
a(3,4)=nu*(1 + nu)*(x1 - x2)*(y1 - y2);
a(4,1)=a(1,4);
a(4,2)=a(2,4);
a(4,3)=a(3,4);
a(4,4)=(1 + nu)*(-(x1 - x2)^2 + (-1 + nu)*(y1 - y2)^2);
a=coef_a*a;

% Local matrix (A sigma, u)^t
b=1/2*[-1 0;0 -1;-1 0;0 -1];

% Local matrix (A sigma, gamma)^t
c=1/4*[y4-y1 x1-x4 y2-y1 x1-x2]';

% Local matrix (A sigma, p)^t
d=zeros(vdim,1);
d(1,1)=(x1 - x4);
d(2,1)=(y1 - y4);
d(3,1)=(x1 - x2);
d(4,1)=(y1 - y2);
d=coef_d*d;

% CC.D en el nodo S-W, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
    % frontera Sur
        Pg1u_L1=int_dir_simpson(i,j,i+1,j,tt,1);
        Pg2u_L1=int_dir_simpson(i,j,i+1,j,tt,2);
    % frontera Oeste
        Pg1u_L4=int_dir_simpson(i,j,i,j+1,tt,1);
        Pg2u_L4=int_dir_simpson(i,j,i,j+1,tt,2);
Pgu=(1/2)*[-Pg1u_L1;-Pg2u_L1;-Pg1u_L4;-Pg2u_L4];  
    
% Non-homogeneous Dir.BC
gamma(pos)=(c'*(a\c))\(c'*(a\Pgu)-(c'*(a\b))*u(ind1u:ind2u)-(c'*(a\d))*p(ind1p));  
pos=pos+1;

% South nodes (j=1)
vdim=6;

% Matrix b is the same for every South node
b=(1/2)*[0 0 -1 0;0 0 0 -1;1 0 -1 0;0 1 0 -1;-1 0 0 0;0 -1 0 0];
a=zeros(vdim,vdim);
d=zeros(vdim,2);

for i=2:N
    ind2u=(i-1+(j-1)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i+(j-1)*N)*2;
%     ind3u=ind4u-1;
    ind1p=ind2u/2;
    ind2p=ind4u/2;
    
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
    denom1=20*JE1r2;
    JE2r1=abs(x5*(-y2 + y3) + x3*(y2 - y5) + x2*(-y3 + y5));
    denom2=20*JE2r1;
    
    a(1,1)=(-((1 + nu)*((-1 + nu)*(x2 - x3)^2 - (y2 - y3)^2)))/denom2;
    a(1,2)=(-(nu*(1 + nu)*(x2 - x3)*(y2 - y3)))/denom2;
    a(1,3)=(-((1 + nu)*((-1 + nu)*(x2 - x3)*(x2 - x5) - (y2 - y3)*(y2 - y5))))/denom2;
    a(1,4)=(-(nu*(1 + nu)*(x2 - x3)*(y2 - y5)))/denom2;
%     a(1,5)=0;
%     a(1,6)=0;
    a(2,1)=a(1,2);
    a(2,2)=(-((1 + nu)*(-(x2 - x3)^2 + (-1 + nu)*(y2 - y3)^2)))/denom2;
    a(2,3)=(-(nu*(1 + nu)*(x2 - x5)*(y2 - y3)))/denom2;
    a(2,4)=(-((1 + nu)*(-((x2 - x3)*(x2 - x5)) + (-1 + nu)*(y2 - y3)*(y2 - y5))))/denom2;
%     a(2,5)=0;
%     a(2,6)=0;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=(-((1 + nu)*((-1 + nu)*(x1 - x2)^2 - (y1 - y2)^2)))/denom1 + ...
           (-((1 + nu)*((-1 + nu)*(x2 - x5)^2 - (y2 - y5)^2)))/denom2;
    a(3,4)=(-(nu*(1 + nu)*(x1 - x2)*(y1 - y2)))/denom1 + ...
           (-(nu*(1 + nu)*(x2 - x5)*(y2 - y5)))/denom2;
    a(3,5)=(-((1 + nu)*((-1 + nu)*(x1 - x2)*(x2 - x3) - (y1 - y2)*(y2 - y3))))/denom1;
    a(3,6)=(-(nu*(1 + nu)*(x1 - x2)*(y2 - y3)))/denom1;
    a(4,1)=a(1,4);
    a(4,2)=a(2,4);
    a(4,3)=a(3,4);
    a(4,4)=((1 + nu)*((x1 - x2)^2 - (-1 + nu)*(y1 - y2)^2))/denom1 + ...
           ((1 + nu)*((x2 - x5)^2 - (-1 + nu)*(y2 - y5)^2))/denom2;
    a(4,5)=(-(nu*(1 + nu)*(x2 - x3)*(y1 - y2)))/denom1;
    a(4,6)=(-((1 + nu)*((-x1 + x2)*(x2 - x3) + (-1 + nu)*(y1 - y2)*(y2 - y3))))/denom1;
%     a(5,1)=0;
%     a(5,2)=0;
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=(-((1 + nu)*((-1 + nu)*(x2 - x3)^2 - (y2 - y3)^2)))/denom1;
    a(5,6)=(-(nu*(1 + nu)*(x2 - x3)*(y2 - y3)))/denom1;
%     a(6,1)=0;
%     a(6,2)=0;
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=(-((1 + nu)*(-(x2 - x3)^2 + (-1 + nu)*(y2 - y3)^2)))/denom1;
    
    c=(1/4)*[y3-y2 x2-x3 y5-y1 x1-x5 y3-y2 x2-x3]';
    
%     d(1,1)=0;
%     d(2,1)=0;
    d(3,1)=(x1 - x2);
    d(4,1)=(y1 - y2);
    d(5,1)=(x2 - x3);
    d(6,1)=(y2 - y3);
    d(1,2)=(x2 - x3);
    d(2,2)=(y2 - y3);
    d(3,2)=(x2 - x5);
    d(4,2)=(y2 - y5);
%     d(5,2)=0;
%     d(6,2)=0;

    d=coef_d*d;
    
    % CC.D en el nodo Sur, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
        % frontera Sur 1
            Pg1u_L2=int_dir_simpson(i,j,i+1,j,tt,1);
            Pg2u_L2=int_dir_simpson(i,j,i+1,j,tt,2);
        % frontera Sur 2
            Pg1u_L1=int_dir_simpson(i-1,j,i,j,tt,1);
            Pg2u_L1=int_dir_simpson(i-1,j,i,j,tt,2);
        Pgu=(1/2)*[-Pg1u_L2;-Pg2u_L2;0;0;-Pg1u_L1;-Pg2u_L1];
    
    % Non-homogeneous Dir.BC
    gamma(pos)=(c'*(a\c))\(c'*(a\Pgu)-(c'*(a\b))*u(ind1u:ind4u)-(c'*(a\d))*p(ind1p:ind2p)); 
    pos=pos+1;
end

% South-East corner node (j=1)
i=N+1;
ind2u=(i-1+(j-1)*N)*2;
ind1u=ind2u-1;
ind1p=ind2u/2;
vdim=4;

x1=x(i-1,j);  
y1=y(i-1,j);            
x2=x(i,j);
y2=y(i,j);
x3=x(i,j+1);
y3=y(i,j+1);
% x4=x(i-1,j+1);
% y4=y(i-1,j+1);

Jer2=abs(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3));
denom=20*Jer2;
coef_a=-1/denom;

a=zeros(vdim,vdim);
a(1,1)=(1 + nu)*((-1 + nu)*(x1 - x2)^2 - (y1 - y2)^2);
a(1,2)=nu*(1 + nu)*(x1 - x2)*(y1 - y2);
a(1,3)=(1 + nu)*((-1 + nu)*(x1 - x2)*(x2 - x3) - (y1 - y2)*(y2 - y3));
a(1,4)=nu*(1 + nu)*(x1 - x2)*(y2 - y3);
a(2,1)=a(1,2);
a(2,2)=(1 + nu)*(-(x1 - x2)^2 + (-1 + nu)*(y1 - y2)^2);
a(2,3)=nu*(1 + nu)*(x2 - x3)*(y1 - y2);
a(2,4)=(1 + nu)*((-x1 + x2)*(x2 - x3) + (-1 + nu)*(y1 - y2)*(y2 - y3));
a(3,1)=a(1,3);
a(3,2)=a(2,3);
a(3,3)=(1 + nu)*((-1 + nu)*(x2 - x3)^2 - (y2 - y3)^2);
a(3,4)=nu*(1 + nu)*(x2 - x3)*(y2 - y3);
a(4,1)=a(1,4);
a(4,2)=a(2,4);
a(4,3)=a(3,4);
a(4,4)=(1 + nu)*(-(x2 - x3)^2 + (-1 + nu)*(y2 - y3)^2);
a=coef_a*a;

b=(1/2)*[1 0;0 1;-1 0;0 -1];

c=(1/4)*[y2-y1 x1-x2 y3-y2 x2-x3]';

d=zeros(vdim,1);
d(1,1)=(x1 - x2);
d(2,1)=(y1 - y2);
d(3,1)=(x2 - x3);
d(4,1)=(y2 - y3);
d=coef_d*d;

% CC.D en el nodo S-E, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
    % frontera Este
        Pg1u_L7=int_dir_simpson(i,j,i,j+1,tt,1);
        Pg2u_L7=int_dir_simpson(i,j,i,j+1,tt,2);
    % frontera Sur
        Pg1u_L3=int_dir_simpson(i-1,j,i,j,tt,1);
        Pg2u_L3=int_dir_simpson(i-1,j,i,j,tt,2);
    Pgu=(1/2)*[Pg1u_L7;Pg2u_L7;-Pg1u_L3;-Pg2u_L3];

% Non-homogeneous Dir.BC
gamma(pos)=(c'*(a\c))\(c'*(a\Pgu)-(c'*(a\b))*u(ind1u:ind2u)-(c'*(a\d))*p(ind1p)); 
pos=pos+1;

for j=2:N

    % West nodes
    i=1;
    ind2u=(i+(j-2)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i+(j-1)*N)*2;
    ind3u=ind4u-1;
    ind1p=ind2u/2;
    ind2p=ind4u/2;
    vdim=6;
    
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
    denom1=20*JE1r4;
    JE2r1=abs(x6*(-y3 + y4) + x4*(y3 - y6) + x3*(-y4 + y6));
    denom2=20*JE2r1;
    
    a=zeros(vdim,vdim);    
    a(1,1)=(-((1 + nu)*((-1 + nu)*(x3 - x4)^2 - (y3 - y4)^2)))/denom1;
    a(1,2)=(-(nu*(1 + nu)*(x3 - x4)*(y3 - y4)))/denom1;
    a(1,3)=((1 + nu)*((-1 + nu)*(x1 - x4)*(x3 - x4) + (y1 - y4)*(-y3 + y4)))/denom1;
    a(1,4)=(nu*(1 + nu)*(x3 - x4)*(y1 - y4))/denom1;
%     a(1,5)=0;
%     a(1,6)=0;
    a(2,1)=a(1,2);
    a(2,2)=(-((1 + nu)*(-(x3 - x4)^2 + (-1 + nu)*(y3 - y4)^2)))/denom1;
    a(2,3)=(nu*(1 + nu)*(x1 - x4)*(y3 - y4))/denom1;
    a(2,4)=((1 + nu)*((x3 - x4)*(-x1 + x4) + (-1 + nu)*(y1 - y4)*(y3 - y4)))/denom1;
%     a(2,5)=0;
%     a(2,6)=0;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=(-((1 + nu)*((-1 + nu)*(x1 - x4)^2 - (y1 - y4)^2)))/denom1 + ...
           (-((1 + nu)*((-1 + nu)*(x4 - x6)^2 - (y4 - y6)^2)))/denom2;
    a(3,4)=(-(nu*(1 + nu)*(x1 - x4)*(y1 - y4)))/denom1 + ...
           (-(nu*(1 + nu)*(x4 - x6)*(y4 - y6)))/denom2;
    a(3,5)=((1 + nu)*((-1 + nu)*(x3 - x4)*(x4 - x6) - (y3 - y4)*(y4 - y6)))/denom2;
    a(3,6)=(nu*(1 + nu)*(x4 - x6)*(y3 - y4))/denom2;
    a(4,1)=a(1,4);
    a(4,2)=a(2,4);
    a(4,3)=a(3,4);
    a(4,4)=((1 + nu)*((x1 - x4)^2 - (-1 + nu)*(y1 - y4)^2))/denom1 + ...
           ((1 + nu)*((x4 - x6)^2 - (-1 + nu)*(y4 - y6)^2))/denom2;
    a(4,5)=(nu*(1 + nu)*(x3 - x4)*(y4 - y6))/denom2;
    a(4,6)=((1 + nu)*((-x3 + x4)*(x4 - x6) + (-1 + nu)*(y3 - y4)*(y4 - y6)))/denom2;
%     a(5,1)=0;
%     a(5,2)=0;
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=(-((1 + nu)*((-1 + nu)*(x3 - x4)^2 - (y3 - y4)^2)))/denom2;
    a(5,6)=(-(nu*(1 + nu)*(x3 - x4)*(y3 - y4)))/denom2;
%     a(6,1)=0;
%     a(6,2)=0;
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=(-((1 + nu)*(-(x3 - x4)^2 + (-1 + nu)*(y3 - y4)^2)))/denom2;
    
    b=(1/2)*[-1 0 0 0;0 -1 0 0;1 0 -1 0;0 1 0 -1;0 0 -1 0;0 0 0 -1];
    
    c=(1/4)*[y3-y4 x4-x3 y6-y1 x1-x6 y3-y4 x4-x3]';
    
    d=zeros(vdim,2);
    d(1,1)=(-x3 + x4);
    d(2,1)=(-y3 + y4);
    d(3,1)=(x1 - x4); 
    d(4,1)=(y1 - y4);
%     d(5,1)=0;
%     d(6,1)=0;
%     d(1,2)=0;
%     d(2,2)=0;
    d(3,2)=(x4 - x6); 
    d(4,2)=(y4 - y6);
    d(5,2)=(x4 - x3);
    d(6,2)=(y4 - y3);
    
    d=coef_d*d;
    
    % CC.D en el nodo Oeste, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
        % frontera Oeste 1
            Pg1u_L4=int_dir_simpson(i,j,i,j-1,tt,1);
            Pg2u_L4=int_dir_simpson(i,j,i,j-1,tt,2);
        % frontera Oeste 2
            Pg1u_L11=int_dir_simpson(i,j,i,j+1,tt,1);
            Pg2u_L11=int_dir_simpson(i,j,i,j+1,tt,2);
        Pgu=(1/2)*[-Pg1u_L4;-Pg2u_L4;0;0;-Pg1u_L11;-Pg2u_L11];
    
    % Non-homogeneous Dir.BC
    gamma(pos)=(c'*(a\c))\(c'*(a\Pgu)-(c'*(a\b))*u([ind1u:ind2u,ind3u:ind4u],1)-(c'*(a\d))*p([ind1p,ind2p],1));
    pos=pos+1;
    
    % Central nodes 
    vdim=8;

    a=zeros(vdim,vdim);
    d=zeros(vdim,4);
    % Matrix b is the same for every central node
    b=(1/2)*[1 0 -1 0 0 0 0 0;0 1 0 -1 0 0 0 0;0 0 1 0 0 0 -1 0;...
             0 0 0 1 0 0 0 -1;0 0 0 0 1 0 -1 0;0 0 0 0 0 1 0 -1;...
             1 0 0 0 -1 0 0 0;0 1 0 0 0 -1 0 0];
         
    for i=2:N
        ind2u=(i-1+(j-2)*N)*2;
        ind1u=ind2u-1;
        ind4u=(i+(j-2)*N)*2;
%         ind3u=ind4u-1;
        ind6u=(i-1+(j-1)*N)*2;
        ind5u=ind6u-1;
        ind8u=(i+(j-1)*N)*2;
%         ind7u=ind8u-1;
        ind1p=ind2u/2;
        ind2p=ind4u/2;
        ind3p=ind6u/2;
        ind4p=ind8u/2;
        
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
%     x7=x(i+1,j+1);
%     y7=y(i+1,j+1);
    x8=x(i,j+1);
    y8=y(i,j+1);
%     x9=x(i-1,j+1);
%     y9=y(i-1,j+1);
    
    JE1r3=abs(x4*(y2 - y3) + x2*(y3 - y4) + x3*(-y2 + y4));
    denom1=20*JE1r3;
    JE2r4=abs(x6*(-y2 + y3) + x3*(y2 - y6) + x2*(-y3 + y6));
    denom2=20*JE2r4;
    JE3r1=abs(x8*(y3 - y6) + x3*(y6 - y8) + x6*(-y3 + y8));
    denom3=20*JE3r1;
    JE4r2=abs(x8*(-y3 + y4) + x4*(y3 - y8) + x3*(-y4 + y8));
    denom4=20*JE4r2;
    
    a(1,1)=(-((1 + nu)*((-1 + nu)*(x3 - x4)^2 - (y3 - y4)^2)))/denom1 + ...
           (-((1 + nu)*((-1 + nu)*(x3 - x6)^2 - (y3 - y6)^2)))/denom2;
    a(1,2)=(-(nu*(1 + nu)*(x3 - x4)*(y3 - y4)))/denom1 + ...
           (-(nu*(1 + nu)*(x3 - x6)*(y3 - y6)))/denom2;
    a(1,3)=(-((1 + nu)*((-1 + nu)*(x2 - x3)*(x3 - x6) - (y2 - y3)*(y3 - y6))))/denom2;
    a(1,4)=(-(nu*(1 + nu)*(x3 - x6)*(y2 - y3)))/denom2;
%     a(1,5)=0;
%     a(1,6)=0;
    a(1,7)=((1 + nu)*((-1 + nu)*(x2 - x3)*(x3 - x4) - (y2 - y3)*(y3 - y4)))/denom1;
    a(1,8)=(nu*(1 + nu)*(x3 - x4)*(y2 - y3))/denom1;
    a(2,1)=a(1,2);
    a(2,2)=((1 + nu)*((x3 - x4)^2 - (-1 + nu)*(y3 - y4)^2))/denom1 + ...
           ((1 + nu)*((x3 - x6)^2 - (-1 + nu)*(y3 - y6)^2))/denom2;
    a(2,3)=(-(nu*(1 + nu)*(x2 - x3)*(y3 - y6)))/denom2;
    a(2,4)=(-((1 + nu)*((-x2 + x3)*(x3 - x6) + (-1 + nu)*(y2 - y3)*(y3 - y6))))/denom2;
%     a(2,5)=0;
%     a(2,6)=0;
    a(2,7)=(nu*(1 + nu)*(x2 - x3)*(y3 - y4))/denom1;
    a(2,8)=((1 + nu)*((-x2 + x3)*(x3 - x4) + (-1 + nu)*(y2 - y3)*(y3 - y4)))/denom1;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=(-((1 + nu)*((-1 + nu)*(x2 - x3)^2 - (y2 - y3)^2)))/denom2 + ...
           (-((1 + nu)*((-1 + nu)*(x3 - x8)^2 - (y3 - y8)^2)))/denom3;
    a(3,4)=(-(nu*(1 + nu)*(x2 - x3)*(y2 - y3)))/denom2 + ...
           (-(nu*(1 + nu)*(x3 - x8)*(y3 - y8)))/denom3;
    a(3,5)=(-((1 + nu)*((-1 + nu)*(x3 - x6)*(x3 - x8) - (y3 - y6)*(y3 - y8))))/denom3;
    a(3,6)=(-(nu*(1 + nu)*(x3 - x8)*(y3 - y6)))/denom3;
%     a(3,7)=0;
%     a(3,8)=0;
    a(4,1)=a(1,4);
    a(4,2)=a(2,4);
    a(4,3)=a(3,4);
    a(4,4)=((1 + nu)*((x2 - x3)^2 - (-1 + nu)*(y2 - y3)^2))/denom2 + ...
           ((1 + nu)*((x3 - x8)^2 - (-1 + nu)*(y3 - y8)^2))/denom3;
    a(4,5)=(-(nu*(1 + nu)*(x3 - x6)*(y3 - y8)))/denom3;
    a(4,6)=((1 + nu)*((x3 - x6)*(x3 - x8) - (-1 + nu)*(y3 - y6)*(y3 - y8)))/denom3;
%     a(4,7)=0;
%     a(4,8)=0;
%     a(5,1)=0;
%     a(5,2)=0;
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=(-((1 + nu)*((-1 + nu)*(x3 - x6)^2 - (y3 - y6)^2)))/denom4 + ...
           (-((1 + nu)*((-1 + nu)*(x3 - x4)^2 - (y3 - y4)^2)))/denom3;
    a(5,6)=(-(nu*(1 + nu)*(x3 - x6)*(y3 - y6)))/denom4 + ....
           (-(nu*(1 + nu)*(x3 - x4)*(y3 - y4)))/denom3;
    a(5,7)=((1 + nu)*((-1 + nu)*(x3 - x4)*(x3 - x8) - (y3 - y4)*(y3 - y8)))/denom4;
    a(5,8)=(nu*(1 + nu)*(x3 - x4)*(y3 - y8))/denom4;
%     a(6,1)=0;
%     a(6,2)=0;
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=((1 + nu)*((x3 - x6)^2 - (-1 + nu)*(y3 - y6)^2))/denom4 + ...
           ((1 + nu)*((x3 - x4)^2 - (-1 + nu)*(y3 - y4)^2))/denom3;
    a(6,7)=(nu*(1 + nu)*(x3 - x8)*(y3 - y4))/denom4;
    a(6,8)=(-((1 + nu)*((x3 - x4)*(x3 - x8) - (-1 + nu)*(y3 - y4)*(y3 - y8))))/denom4;
    a(7,1)=a(1,7);
    a(7,2)=a(2,7);
%     a(7,3)=0;
%     a(7,4)=0;
    a(7,5)=a(5,7);
    a(7,6)=a(6,7);
    a(7,7)=(-((1 + nu)*((-1 + nu)*(x2 - x3)^2 - (y2 - y3)^2)))/denom1 + ...
           (-((1 + nu)*((-1 + nu)*(x3 - x8)^2 - (y3 - y8)^2)))/denom4;
    a(7,8)=(-(nu*(1 + nu)*(x2 - x3)*(y2 - y3)))/denom1 + ... 
           (-(nu*(1 + nu)*(x3 - x8)*(y3 - y8)))/denom4;
    a(8,1)=a(1,8);
    a(8,2)=a(2,8);
%     a(8,3)=0;
%     a(8,4)=0;
    a(8,5)=a(5,8);
    a(8,6)=a(6,8);
    a(8,7)=a(7,8);
    a(8,8)=((1 + nu)*((x2 - x3)^2 - (-1 + nu)*(y2 - y3)^2))/denom1 + ...
           ((1 + nu)*((x3 - x8)^2 - (-1 + nu)*(y3 - y8)^2))/denom4;
    
    c=(1/4)*[y6-y4 x4-x6 y8-y2 x2-x8 y6-y4 x4-x6 y8-y2 x2-x8]';
    
    % ¡¡OJO con las posiciones 3 y 4 (van intercaladas <- va antes la 4 que 
    % la 3!!
    d(1,1)=(-x3 + x4);
    d(2,1)=(-y3 + y4);
%     d(3,1)=0;
%     d(4,1)=0;
%     d(5,1)=0;
%     d(6,1)=0;
    d(7,1)=(x2 - x3);
    d(8,1)=(y2 - y3);
    d(1,2)=(-x6 + x3);
    d(2,2)=(-y6 + y3);
    d(3,2)=(x2 - x3);
    d(4,2)=(y2 - y3);
%     d(5,2)=0;
%     d(6,2)=0;
%     d(7,2)=0;
%     d(8,2)=0;

%     ¡OJO!: ¡Posicones 3 y 4 van intercaladas! (P4 va antes que P3)

%     d(1,3)=0;
%     d(2,3)=0;
%     d(3,3)=0;
%     d(4,3)=0;
    d(5,3)=(x4 - x3);
    d(6,3)=(y4 - y3);
    d(7,3)=(x3 - x8);
    d(8,3)=(y3 - y8);

%     d(1,4)=0;
%     d(2,4)=0;
    d(3,4)=(x3 - x8);
    d(4,4)=(y3 - y8);
    d(5,4)=(x3 - x6);
    d(6,4)=(y3 - y6);
%     d(7,4)=0;
%     d(8,4)=0;
    
    d=coef_d*d;
    
    % NO Dir.BC
    gamma(pos)=(c'*(a\c))\(-(c'*(a\b))*u([ind1u:ind4u,ind5u:ind8u],1)-(c'*(a\d))*p([ind1p:ind2p,ind3p:ind4p],1)); 
    pos=pos+1;
    end
    
    % East nodes
    i=N+1;
    ind2u=(i-1+(j-2)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i-1+(j-1)*N)*2;
    ind3u=ind4u-1;
    ind1p=ind2u/2;
    ind2p=ind4u/2;

    vdim=6;
    
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
    denom1=20*JE1r3;
    JE2r2=abs(x5*(-y3 + y4) + x4*(y3 - y5) + x3*(-y4 + y5));
    denom2=20*JE2r2;
    
    a=zeros(vdim,vdim);
    a(1,1)=(-((1 + nu)*((-1 + nu)*(x3 - x4)^2 - (y3 - y4)^2)))/denom1;
    a(1,2)=(-(nu*(1 + nu)*(x3 - x4)*(y3 - y4)))/denom1;
%     a(1,3)=0;
%     a(1,4)=0;
    a(1,5)=((1 + nu)*((-1 + nu)*(x2 - x3)*(x3 - x4) - (y2 - y3)*(y3 - y4)))/denom1;
    a(1,6)=(nu*(1 + nu)*(x3 - x4)*(y2 - y3))/denom1;
    a(2,1)=a(1,2);
    a(2,2)=(-((1 + nu)*(-(x3 - x4)^2 + (-1 + nu)*(y3 - y4)^2)))/denom1;
%     a(2,3)=0;
%     a(2,4)=0:
    a(2,5)=(nu*(1 + nu)*(x2 - x3)*(y3 - y4))/denom1;
    a(2,6)=((1 + nu)*((-x2 + x3)*(x3 - x4) + (-1 + nu)*(y2 - y3)*(y3 - y4)))/denom1;
%     a(3,1)=0;
%     a(3,2)=0;
    a(3,3)=(-((1 + nu)*((-1 + nu)*(x3 - x4)^2 - (y3 - y4)^2)))/denom2;
    a(3,4)=(-(nu*(1 + nu)*(x3 - x4)*(y3 - y4)))/denom2;
    a(3,5)=((1 + nu)*((-1 + nu)*(x3 - x4)*(x3 - x5) - (y3 - y4)*(y3 - y5)))/denom2;
    a(3,6)=(nu*(1 + nu)*(x3 - x4)*(y3 - y5))/denom2;
%     a(4,1)=0;
%     a(4,2)=0;
    a(4,3)=a(3,4);
    a(4,4)=(-((1 + nu)*(-(x3 - x4)^2 + (-1 + nu)*(y3 - y4)^2)))/denom2;
    a(4,5)=(nu*(1 + nu)*(x3 - x5)*(y3 - y4))/denom2;
    a(4,6)=((1 + nu)*(-((x3 - x4)*(x3 - x5)) + (-1 + nu)*(y3 - y4)*(y3 - y5)))/denom2;
    a(5,1)=a(1,5);
    a(5,2)=a(2,5);
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=(-((1 + nu)*((-1 + nu)*(x2 - x3)^2 - (y2 - y3)^2)))/denom1 + ...
           (-((1 + nu)*((-1 + nu)*(x3 - x5)^2 - (y3 - y5)^2)))/denom2;
    a(5,6)=(-(nu*(1 + nu)*(x2 - x3)*(y2 - y3)))/denom1 + ...
           (-(nu*(1 + nu)*(x3 - x5)*(y3 - y5)))/denom2;
    a(6,1)=a(1,6);
    a(6,2)=a(2,6);
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=((1 + nu)*((x2 - x3)^2 - (-1 + nu)*(y2 - y3)^2))/denom1 +...
           ((1 + nu)*((x3 - x5)^2 - (-1 + nu)*(y3 - y5)^2))/denom2;
    
    b=(1/2)*[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1;1 0 -1 0;0 1 0 -1];
    
    c=(1/4)*[y3-y4 x4-x3 y3-y4 x4-x3 y5-y2 x2-x5]';
    
    d=zeros(vdim,2);
    d(1,1)=(-x3 + x4);
    d(2,1)=(-y3 + y4);
%     d(3,1)=0;
%     d(4,1)=0;
    d(5,1)=(x2 - x3);
    d(6,1)=(y2 - y3);
%     d(1,2)=0;
%     d(2,2)=0;
    d(3,2)=(x4 - x3);
    d(4,2)=(y4 - y3);
    d(5,2)=(x3 - x5);
    d(6,2)=(y3 - y5);
    
    d=coef_d*d;
    
    % CC.D en el nodo Este, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
        % frontera Este 1
            Pg1u_L7=int_dir_simpson(i,j,i,j-1,tt,1);
            Pg2u_L7=int_dir_simpson(i,j,i,j-1,tt,2);
        % frontera Este 2
            Pg1u_L14=int_dir_simpson(i,j,i,j+1,tt,1);
            Pg2u_L14=int_dir_simpson(i,j,i,j+1,tt,2);
        Pgu=(1/2)*[Pg1u_L7;Pg2u_L7;Pg1u_L14;Pg2u_L14;0;0];
    
    % Non-homogeneous Dir.BC
    gamma(pos)=(c'*(a\c))\(c'*(a\Pgu)-(c'*(a\b))*u([ind1u:ind2u,ind3u:ind4u],1)-(c'*(a\d))*p([ind1p,ind2p],1)); 
    pos=pos+1;
end

% North-West corner node
i=1;
j=N+1;
ind2u=(i+(j-2)*N)*2;
ind1u=ind2u-1;
ind1p=ind2u/2;
vdim=4;

x1=x(i,j-1);
y1=y(i,j-1);           
% x2=x(i+1,j-1);
% y2=y(i+1,j-1);
x3=x(i+1,j);
y3=y(i+1,j);
x4=x(i,j);
y4=y(i,j);

Jer4=abs(x4*(y1 - y3) + x1*(y3 - y4) + x3*(-y1 + y4));
denom=20*Jer4;

a=zeros(vdim,vdim);
a(1,1)=-((1 + nu)*((-1 + nu)*(x3 - x4)^2 - (y3 - y4)^2));
a(1,2)=-(nu*(1 + nu)*(x3 - x4)*(y3 - y4));
a(1,3)=(1 + nu)*((-1 + nu)*(x1 - x4)*(x3 - x4) + (y1 - y4)*(-y3 + y4));
a(1,4)=nu*(1 + nu)*(x3 - x4)*(y1 - y4);
a(2,1)=a(1,2);
a(2,2)=-((1 + nu)*(-(x3 - x4)^2 + (-1 + nu)*(y3 - y4)^2));
a(2,3)=nu*(1 + nu)*(x1 - x4)*(y3 - y4);
a(2,4)=(1 + nu)*((x3 - x4)*(-x1 + x4) + (-1 + nu)*(y1 - y4)*(y3 - y4));
a(3,1)=a(1,3);
a(3,2)=a(2,3);
a(3,3)=-((1 + nu)*((-1 + nu)*(x1 - x4)^2 - (y1 - y4)^2));
a(3,4)=-(nu*(1 + nu)*(x1 - x4)*(y1 - y4));
a(4,1)=a(1,4);
a(4,2)=a(2,4);
a(4,3)=a(3,4);
a(4,4)=-((1 + nu)*(-(x1 - x4)^2 + (-1 + nu)*(y1 - y4)^2));
a=(1/denom)*a;

b=(1/2)*[-1 0;0 -1;1 0;0 1];

c=(1/4)*[y3-y4 x4-x3 y4-y1 x1-x4]';

d=zeros(vdim,1);
d(1,1)=(-x3 + x4);
d(2,1)=(-y3 + y4);
d(3,1)=(x1 - x4);
d(4,1)=(y1 - y4);
d=coef_d*d;

% CC.D en el nodo N-W, fórmula de cuadratura para los término Pg1 y Pg2 en:
    % frontera Oeste
        Pg1u_L18=int_dir_simpson(i,j,i,j-1,tt,1);
        Pg2u_L18=int_dir_simpson(i,j,i,j-1,tt,2);
    % frontera Norte
        Pg1u_L22=int_dir_simpson(i,j,i+1,j,tt,1);
        Pg2u_L22=int_dir_simpson(i,j,i+1,j,tt,2);
    Pgu=(1/2)*[-Pg1u_L18;-Pg2u_L18;Pg1u_L22;Pg2u_L22];
    
% Non-homogeneous Dir.BC
gamma(pos)=(c'*(a\c))\(c'*(a\Pgu)-(c'*(a\b))*u(ind1u:ind2u)-(c'*(a\d))*p(ind1p)); 
pos=pos+1;

% North nodes (j=N+1)
vdim=6;

% Matrix b is the same for every North node
b=(1/2)*[1 0 -1 0;0 1 0 -1;0 0 1 0;0 0 0 1;1 0 0 0;0 1 0 0];

a=zeros(vdim,vdim);
d=zeros(vdim,2);

for i=2:N
    ind2u=(i-1+(j-2)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i+(j-2)*N)*2;
%     ind3u=ind4u-1;
    ind1p=ind2u/2;
    ind2p=ind4u/2;
    
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
    denom1=20*JE1r3;
    JE2r4=abs(x6*(-y2 + y3) + x3*(y2 - y6) + x2*(-y3 + y6));
    denom2=20*JE2r4;
    
    a(1,1)=(-((1 + nu)*((-1 + nu)*(x3 - x4)^2 - (y3 - y4)^2)))/denom1 + ...
           (-((1 + nu)*((-1 + nu)*(x3 - x6)^2 - (y3 - y6)^2)))/denom2;
    a(1,2)=(-(nu*(1 + nu)*(x3 - x4)*(y3 - y4)))/denom1 + ...
           (-(nu*(1 + nu)*(x3 - x6)*(y3 - y6)))/denom2;
    a(1,3)=(-((1 + nu)*((-1 + nu)*(x2 - x3)*(x3 - x6) - (y2 - y3)*(y3 - y6))))/denom2;
    a(1,4)=(-(nu*(1 + nu)*(x3 - x6)*(y2 - y3)))/denom2;
    a(1,5)=((1 + nu)*((-1 + nu)*(x2 - x3)*(x3 - x4) - (y2 - y3)*(y3 - y4)))/denom1;
    a(1,6)=(nu*(1 + nu)*(x3 - x4)*(y2 - y3))/denom1;
    a(2,1)=a(1,2);
    a(2,2)=((1 + nu)*((x3 - x4)^2 - (-1 + nu)*(y3 - y4)^2))/denom1 + ...
           ((1 + nu)*((x3 - x6)^2 - (-1 + nu)*(y3 - y6)^2))/denom2;
    a(2,3)=(-(nu*(1 + nu)*(x2 - x3)*(y3 - y6)))/denom2;
    a(2,4)=(-((1 + nu)*((-x2 + x3)*(x3 - x6) + (-1 + nu)*(y2 - y3)*(y3 - y6))))/denom2;
    a(2,5)=(nu*(1 + nu)*(x2 - x3)*(y3 - y4))/denom1;
    a(2,6)=((1 + nu)*((-x2 + x3)*(x3 - x4) + (-1 + nu)*(y2 - y3)*(y3 - y4)))/denom1;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=(-((1 + nu)*((-1 + nu)*(x2 - x3)^2 - (y2 - y3)^2)))/denom2;
    a(3,4)=(-(nu*(1 + nu)*(x2 - x3)*(y2 - y3)))/denom2;
%     a(3,5)=0;
%     a(3,6)=0;
    a(4,1)= a(1,4);
    a(4,2)= a(2,4);
    a(4,3)= a(3,4);
    a(4,4)=(-((1 + nu)*(-(x2 - x3)^2 + (-1 + nu)*(y2 - y3)^2)))/denom2;
%     a(4,5)=0;
%     a(4,6)=0;
    a(5,1)=a(1,5);
    a(5,2)=a(2,5);
%     a(5,3)=0;
%     a(5,4)=0;
    a(5,5)=(-((1 + nu)*((-1 + nu)*(x2 - x3)^2 - (y2 - y3)^2)))/denom1;
    a(5,6)=(-(nu*(1 + nu)*(x2 - x3)*(y2 - y3)))/denom1;
    a(6,1)=a(1,6);
    a(6,2)=a(2,6);
%     a(6,3)=0;
%     a(6,4)=0;
    a(6,5)=a(5,6);
    a(6,6)=(-((1 + nu)*(-(x2 - x3)^2 + (-1 + nu)*(y2 - y3)^2)))/denom1;
    
    c=(1/4)*[y6-y4 x4-x6 y3-y2 x2-x3 y3-y2 x2-x3]';
    
    d(1,1)=(-x3 + x4);
    d(2,1)=(-y3 + y4);
%     d(3,1)=0;
%     d(4,1)=0;
    d(5,1)=(x2 - x3);
    d(6,1)=(y2 - y3);
    d(1,2)=(-x6 + x3);
    d(2,2)=(-y6 + y3);
    d(3,2)=(x2 - x3);
    d(4,2)=(y2 - y3);
%     d(5,2)=0;
%     d(6,2)=0;

    d=coef_d*d;
    
    % CC.D en el nodo Norte, fórmula de cuadratura para los términos Pg1 y Pg2 de u en:
        % frontera Norte 1
            Pg1u_L23=int_dir_simpson(i,j,i+1,j,tt,1);
            Pg2u_L23=int_dir_simpson(i,j,i+1,j,tt,2);
        % frontera Norte 2
            Pg1u_L22=int_dir_simpson(i-1,j,i,j,tt,1);
            Pg2u_L22=int_dir_simpson(i-1,j,i,j,tt,2);
        Pgu=(1/2)*[0;0;Pg1u_L23;Pg2u_L23;Pg1u_L22;Pg2u_L22];
        
    % Non-homogeneous Dir.BC
    gamma(pos)=(c'*(a\c))\(c'*(a\Pgu)-(c'*(a\b))*u(ind1u:ind4u)-(c'*(a\d))*p(ind1p:ind2p));
    pos=pos+1;
end

% North-East corner node (j=N+1)
i=N+1;
ind2u=(i-1+(j-2)*N)*2;
ind1u=ind2u-1;
ind1p=ind2u/2;
vdim=4;

% x1=x(i-1,j-1);  
% y1=y(i-1,j-1);            
x2=x(i,j-1);
y2=y(i,j-1);
x3=x(i,j);
y3=y(i,j);
x4=x(i-1,j);
y4=y(i-1,j);

Jer3=abs(x4*(y2 - y3) + x2*(y3 - y4) + x3*(-y2 + y4));
denom=20*Jer3;

a=zeros(vdim,vdim);
a(1,1)=-((1 + nu)*((-1 + nu)*(x3 - x4)^2 - (y3 - y4)^2));
a(1,2)=-(nu*(1 + nu)*(x3 - x4)*(y3 - y4));
a(1,3)=(1 + nu)*((-1 + nu)*(x2 - x3)*(x3 - x4) - (y2 - y3)*(y3 - y4));
a(1,4)=nu*(1 + nu)*(x3 - x4)*(y2 - y3);
a(2,1)=a(1,2);
a(2,2)=-((1 + nu)*(-(x3 - x4)^2 + (-1 + nu)*(y3 - y4)^2));
a(2,3)=nu*(1 + nu)*(x2 - x3)*(y3 - y4);
a(2,4)=(1 + nu)*((-x2 + x3)*(x3 - x4) + (-1 + nu)*(y2 - y3)*(y3 - y4));
a(3,1)=a(1,3);
a(3,2)=a(2,3);
a(3,3)=-((1 + nu)*((-1 + nu)*(x2 - x3)^2 - (y2 - y3)^2));
a(3,4)=-(nu*(1 + nu)*(x2 - x3)*(y2 - y3));
a(4,1)=a(1,4);
a(4,2)=a(2,4);
a(4,3)=a(3,4);
a(4,4)=-((1 + nu)*(-(x2 - x3)^2 + (-1 + nu)*(y2 - y3)^2));
a=(1/denom)*a;

b=(1/2)*[1 0;0 1;1 0;0 1];

c=(1/4)*[y3-y4 x4-x3 y3-y2 x2-x3]';

d=zeros(vdim,1);
d(1,1)=(-x3 + x4);
d(2,1)=(-y3 + y4);
d(3,1)=(x2 - x3);
d(4,1)=(y2 - y3);
d=coef_d*d;

% CC.D en el nodo N-E, fórmula de cuadratura para los términos Pg1 y Pg2 de u en:
    % frontera Este
        Pg1u_L21=int_dir_simpson(i,j,i,j-1,tt,1);
        Pg2u_L21=int_dir_simpson(i,j,i,j-1,tt,2);
    % frontera Norte
        Pg1u_L24=int_dir_simpson(i-1,j,i,j,tt,1);
        Pg2u_L24=int_dir_simpson(i-1,j,i,j,tt,2);
    Pgu=(1/2)*[Pg1u_L21;Pg2u_L21;Pg1u_L24;Pg2u_L24];

% Non-homogeneous Dir.BC
gamma(pos)=(c'*(a\c))\(c'*(a\Pgu)-(c'*(a\b))*u(ind1u:ind2u)-(c'*(a\d))*p(ind1p)); 
return
end
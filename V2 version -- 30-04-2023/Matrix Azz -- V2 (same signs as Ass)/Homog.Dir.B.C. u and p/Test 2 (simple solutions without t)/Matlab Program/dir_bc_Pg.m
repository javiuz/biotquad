function [gDu,gDp]=dir_bc_Pg(delta_t,tt)

global NN x y alpha lambda mu Kinv

N=NN;

gDu=zeros(2*N*N,1);
gDp=zeros(N*N,1);

coef_denom_elas=(16.*mu*(lambda + mu));
coef_d=alpha/(8*(lambda+mu));

% South-West corner node
i=1;
j=1;
ind2u=(i+(j-1)*N)*2;
ind1u=ind2u-1;
ind1p=ind2u/2;
vdim=4;
vdim2=vdim/2;

x1=x(i,j);  
y1=y(i,j);            
x2=x(i+1,j);
y2=y(i+1,j);
% x3=x(i+1,j+1);
% y3=y(i+1,j+1);
x4=x(i,j+1);
y4=y(i,j+1);

Jer1=abs(x2*y1 - x4*y1 - x1*y2 + x4*y2 + x1*y4 - x2*y4);
SW_denom_ss=Jer1*coef_denom_elas;
SW_denom_zz=Jer1*4;

% Matriz local (A sigma, sigma)
a=zeros(vdim,vdim);
a(1,1)=((lambda + 2*mu)*(x1 - x4)^2 + 2*(lambda + mu)*(y1 - y4)^2);
a(2,1)=(lambda*(-x1 + x4)*(y1 - y4));
a(3,1)=((lambda + 2*mu)*(x1 - x2)*(x1 - x4) + 2*(lambda + mu)*(y1 - y2)*(y1 - y4));
a(4,1)=-(lambda*(x1 - x4)*(y1 - y2));
a(1,2)=a(2,1);
a(2,2)=(2*(lambda + mu)*(x1 - x4)^2 + (lambda + 2*mu)*(y1 - y4)^2);
a(3,2)=-(lambda*(x1 - x2)*(y1 - y4));
a(4,2)=(2*(lambda + mu)*(x1 - x2)*(x1 - x4) + (lambda + 2*mu)*(y1 - y2)*(y1 - y4));
a(1,3)=a(3,1);
a(2,3)=a(3,2);
a(3,3)=((lambda + 2*mu)*(x1 - x2)^2 + 2*(lambda + mu)*(y1 - y2)^2);
a(4,3)=-(lambda*(x1 - x2)*(y1 - y2));
a(1,4)=a(4,1);
a(2,4)=a(4,2);
a(3,4)=a(4,3);
a(4,4)=(2*(lambda + mu)*(x1 - x2)^2 + (lambda + 2*mu)*(y1 - y2)^2);
a=a/SW_denom_ss;

% Matriz local (A sigma, u)^t
b=zeros(vdim,2);
b(1,1)=-(1/2);
% b(1,2)=0;
% b(2,1)=0;
b(2,2)=-(1/2);
b(3,1)=-(1/2);
% b(3,2)=0;
% b(4,1)=0;
b(4,2)=-(1/2);

% Matriz local (A sigma, gamma)^t
c=(1/4)*[(-y1 + y4),(x1 - x4),(-y1 + y2),(x1 - x2)]';

% Matriz local (A sigma, p)^t
d=zeros(vdim,1);
d(1,1)=(-x1 + x4);
d(2,1)=(-y1 + y4);
d(3,1)=(-x1 + x2);
d(4,1)=(-y1 + y2);
d=coef_d*d;

% Matriz local (A z, p)^T
f=zeros(vdim2,1);
f(1,1)=1/2;
f(2,1)=1/2;

% Matriz local (A z, z)^{-1}
t=zeros(vdim2,vdim2);
t(1,1)=Kinv(1,1)*(x1 - x4)^2 + (2*Kinv(1,2)*(x1 - x4) + Kinv(2,2)*(y1 - y4))*(y1 - y4);
t(2,1)=(x1 - x4)*(Kinv(1,1)*(x1 - x2) + Kinv(1,2)*(y1 - y2)) + ...
       (Kinv(1,2)*(x1 - x2) + Kinv(2,2)*(y1 - y2))*(y1 - y4);
t(1,2)=t(2,1);
t(2,2)=Kinv(1,1)*(x1 - x2)^2 + (2*Kinv(1,2)*(x1 - x2) + Kinv(2,2)*(y1 - y2))*(y1 - y2);
t=t/SW_denom_zz;

% CC.D en el nodo S-W, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
    % frontera Sur
        Pg1u_L1=int_dir_simpson(i,j,i+1,j,tt,1);
        Pg2u_L1=int_dir_simpson(i,j,i+1,j,tt,2);
    % frontera Oeste
        Pg1u_L4=int_dir_simpson(i,j,i,j+1,tt,1);
        Pg2u_L4=int_dir_simpson(i,j,i,j+1,tt,2);
        
    Pgu=(1/2)*[-Pg1u_L1;-Pg2u_L1;-Pg1u_L4;-Pg2u_L4];  
    localv_u=(b'-((b'*(a\c))*((c'*(a\c))\c')))*(a\Pgu);
    gDu(ind1u:ind2u)=gDu(ind1u:ind2u)+localv_u;
   
% CC.D en el nodo S-W, fórmula de cuadratura para el término Pg de p en:
    % frontera Sur
        Pgp_L1=int_dir_simpson(i,j,i+1,j,tt,3);
    % frontera Oeste
        Pgp_L4=int_dir_simpson(i,j,i,j+1,tt,3);
    
        Pgp=(1/2)*[Pgp_L1;Pgp_L4];
        localv_p=((delta_t*f')*(t\Pgp))+((-d'+((d'*(a\c))*((c'*(a\c))\c')))*(a\Pgu));
        gDp(ind1p)=gDp(ind1p)+localv_p;

% South nodes (j=1)
vdim=6;
vdim2=vdim/2;

% Matrices b and f are the same for every South node
b=zeros(vdim,4);
b(1,3)=-(1/2);
% b(1,4)=0;
% b(2,3)=0;
b(2,4)=-(1/2);
b(3,1)=1/2;
% b(3,2)=0;
b(3,3)=-(1/2);
% b(3,4)=0;
% b(4,1)=0;
b(4,2)=1/2;
% b(4,3)=0;
b(4,4)=-(1/2);
b(5,1)=-(1/2);
% b(5,2)=0;
% b(6,1)=0;
b(6,2)=-(1/2);

f=zeros(vdim2,2);
% f(1,1)=0;
f(2,1)=-(1/2);
f(3,1)=1/2;
f(1,2)=1/2;
f(2,2)=1/2;
% f(3,2)=0;

a=zeros(vdim,vdim);
d=zeros(vdim,2);
t=zeros(vdim2,vdim2);

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
    
    JE1r2=abs(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3);
    S_denom_E1_ss=coef_denom_elas*JE1r2;
    S_denom_E1_zz=4*JE1r2;
    JE2r1=abs(x5*(-y2 + y3) + x3*(y2 - y5) + x2*(-y3 + y5));
    S_denom_E2_ss=coef_denom_elas*JE2r1;
    S_denom_E2_zz=4*JE2r1;
    
    a(1,1)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/S_denom_E2_ss;
    a(2,1)=-(lambda*(x2 - x3)*(y2 - y3))/S_denom_E2_ss;
    a(3,1)=((lambda + 2*mu)*(x2 - x3)*(x2 - x5) + 2*(lambda + mu)*(y2 - y3)*(y2 - y5))/S_denom_E2_ss;
    a(4,1)=-(lambda*(x2 - x3)*(y2 - y5))/S_denom_E2_ss;
%     a(5,1)=0;
%     a(6,1)=0;
    a(1,2)=a(2,1);
    a(2,2)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/S_denom_E2_ss;
    a(3,2)=-(lambda*(x2 - x5)*(y2 - y3))/S_denom_E2_ss;
    a(4,2)=(2*(lambda + mu)*(x2 - x3)*(x2 - x5) + (lambda + 2*mu)*(y2 - y3)*(y2 - y5))/S_denom_E2_ss;
%     a(5,2)=0;
%     a(6,2)=0;
    a(1,3)=a(3,1);
    a(2,3)=a(3,2);
    a(3,3)=((lambda + 2*mu)*(x1 - x2)^2 + 2*(lambda + mu)*(y1 - y2)^2)/S_denom_E1_ss+...
           ((lambda + 2*mu)*(x2 - x5)^2 + 2*(lambda + mu)*(y2 - y5)^2)/S_denom_E2_ss;
    a(4,3)=-(lambda*(x1 - x2)*(y1 - y2))/S_denom_E1_ss - ...
           (lambda*(x2 - x5)*(y2 - y5))/S_denom_E2_ss;
    a(5,3)=((lambda + 2*mu)*(x1 - x2)*(x2 - x3) + 2*(lambda + mu)*(y1 - y2)*(y2 - y3))/S_denom_E1_ss;
    a(6,3)=-(lambda*(x1 - x2)*(y2 - y3))/S_denom_E1_ss;
    a(1,4)=a(4,1);
    a(2,4)=a(4,2);
    a(3,4)=a(4,3);
    a(4,4)=(2*(lambda + mu)*(x1 - x2)^2 + (lambda + 2*mu)*(y1 - y2)^2)/S_denom_E1_ss+...
           (2*(lambda + mu)*(x2 - x5)^2 + (lambda + 2*mu)*(y2 - y5)^2)/S_denom_E2_ss;
    a(5,4)=-(lambda*(x2 - x3)*(y1 - y2))/S_denom_E1_ss;
    a(6,4)=(2*(lambda + mu)*(x1 - x2)*(x2 - x3) + (lambda + 2*mu)*(y1 - y2)*(y2 - y3))/S_denom_E1_ss;
%     a(1,5)=0;
%     a(2,5)=0;
    a(3,5)=a(5,3);
    a(4,5)=a(5,4);
    a(5,5)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/S_denom_E1_ss;
    a(6,5)=-(lambda*(x2 - x3)*(y2 - y3))/S_denom_E1_ss;
%     a(1,6)=0;
%     a(2,6)=0;
    a(3,6)=a(6,3);
    a(4,6)=a(6,4);
    a(5,6)=a(6,5);
    a(6,6)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/S_denom_E1_ss;
    
    c=(1/4)*[(-y2 + y3),(x2 - x3),(-y1 + y5),(x1 - x5),(-y2 + y3),(x2 - x3)]';
    
%     d(1,1)=0;
%     d(2,1)=0;
    d(3,1)=(-x1 + x2);
    d(4,1)=(-y1 + y2);
    d(5,1)=(-x2 + x3);
    d(6,1)=(-y2 + y3);
    d(1,2)=(-x2 + x3);
    d(2,2)=(-y2 + y3);
    d(3,2)=(-x2 + x5);
    d(4,2)=(-y2 + y5);
%     d(5,2)=0;
%     d(6,2)=0;    
    d=coef_d*d;

    t(1,1)=(Kinv(1,1)*(x2 - x3)^2 + (2*Kinv(1,2)*(x2 - x3) + Kinv(2,2)*(y2 - y3))*(y2 - y3))/S_denom_E2_zz;
    t(2,1)=((y2 - y3)*(Kinv(1,2)*(x2 - x5) + Kinv(2,2)*(y2 - y5)) + ...
            (x2 - x3)*(Kinv(1,1)*(x2 - x5) + Kinv(1,2)*(y2 - y5)))/S_denom_E2_zz;
%     t(3,1)=0;
    t(1,2)=t(2,1);
    t(2,2)=(Kinv(1,1)*(x1 - x2)^2 + (2*Kinv(1,2)*(x1 - x2) + Kinv(2,2)*(y1 - y2))*(y1 - y2))/S_denom_E1_zz +...
           (Kinv(1,1)*(x2 - x5)^2 + (2*Kinv(1,2)*(x2 - x5) + Kinv(2,2)*(y2 - y5))*(y2 - y5))/S_denom_E2_zz;
    t(3,2)=(Kinv(1,1)*(x1 - x2)*(x2 - x3) + Kinv(2,2)*(y1 - y2)*(y2 - y3) + Kinv(1,2)*(x3*(-y1 + y2) + x1*(y2 - y3) + x2*(y1 - 2*y2 + y3)))/S_denom_E1_zz;
%     t(1,3)=0;
    t(2,3)=t(3,2);
    t(3,3)=(Kinv(1,1)*(x2 - x3)^2 + (2*Kinv(1,2)*(x2 - x3) + Kinv(2,2)*(y2 - y3))*(y2 - y3))/S_denom_E1_zz;
    
    % CC.D en el nodo Sur, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
        % frontera Sur 1
            Pg1u_L2=int_dir_simpson(i,j,i+1,j,tt,1);
            Pg2u_L2=int_dir_simpson(i,j,i+1,j,tt,2);
        % frontera Sur 2
            Pg1u_L1=int_dir_simpson(i-1,j,i,j,tt,1);
            Pg2u_L1=int_dir_simpson(i-1,j,i,j,tt,2);
        Pgu=(1/2)*[-Pg1u_L2;-Pg2u_L2;0;0;-Pg1u_L1;-Pg2u_L1];
    
        localv_u=(b'-((b'*(a\c))*((c'*(a\c))\c')))*(a\Pgu);
        gDu(ind1u:ind4u)=gDu(ind1u:ind4u)+localv_u;
        
    % CC.D en los nodos Sur, fórmula de cuadratura para el término Pg de p en:
        % frontera Sur 1
            Pgp_L2=int_dir_simpson(i,j,i+1,j,tt,3); 
        % frontera Sur 2
            Pgp_L1=int_dir_simpson(i-1,j,i,j,tt,3); 
        Pgp=(1/2)*[Pgp_L2;0;Pgp_L1];
        
        localv_p=((delta_t*f')*(t\Pgp))+((-d'+((d'*(a\c))*((c'*(a\c))\c')))*(a\Pgu));
        gDp(ind1p:ind2p)=gDp(ind1p:ind2p)+localv_p;
end

% South-East corner node (j=1)
i=N+1;
ind2u=(i-1+(j-1)*N)*2;
ind1u=ind2u-1;
ind1p=ind2u/2;
vdim=4;
vdim2=vdim/2;

x1=x(i-1,j);  
y1=y(i-1,j);            
x2=x(i,j);
y2=y(i,j);
x3=x(i,j+1);
y3=y(i,j+1);
% x4=x(i-1,j+1);
% y4=y(i-1,j+1);

Jer2=abs(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3);
SE_denom_ss=Jer2*coef_denom_elas;
SE_denom_zz=Jer2*4;

a=zeros(vdim,vdim);
a(1,1)=((lambda + 2*mu)*(x1 - x2)^2 + 2*(lambda + mu)*(y1 - y2)^2);
a(2,1)=-(lambda*(x1 - x2)*(y1 - y2));
a(3,1)=((lambda + 2*mu)*(x1 - x2)*(x2 - x3) + 2*(lambda + mu)*(y1 - y2)*(y2 - y3));
a(4,1)=-(lambda*(x1 - x2)*(y2 - y3));
a(1,2)=a(2,1);
a(2,2)=(2*(lambda + mu)*(x1 - x2)^2 + (lambda + 2*mu)*(y1 - y2)^2);
a(3,2)=-(lambda*(x2 - x3)*(y1 - y2));
a(4,2)=(2*(lambda + mu)*(x1 - x2)*(x2 - x3) + (lambda + 2*mu)*(y1 - y2)*(y2 - y3));
a(1,3)=a(3,1);
a(2,3)=a(3,2);
a(3,3)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2);
a(4,3)=-(lambda*(x2 - x3)*(y2 - y3));
a(1,4)=a(4,1);
a(2,4)=a(4,2);
a(3,4)=a(4,3);
a(4,4)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2);
a=a/SE_denom_ss;

b=zeros(vdim,2);
b(1,1)=1/2;
% b(1,2)=0;
% b(2,1)=0;
b(2,2)=1/2;
b(3,1)=-(1/2);
% b(3,2)=0;
% b(4,1)=0;
b(4,2)=-(1/2);

c=(1/4)*[(-y1 + y2),(x1 - x2),(-y2 + y3),(x2 - x3)]';

d=zeros(vdim,1);
d(1,1)=(-x1 + x2);
d(2,1)=(-y1 + y2);
d(3,1)=(-x2 + x3);
d(4,1)=(-y2 + y3);
d=coef_d*d;

f=zeros(vdim2,1);
f(1,1)=-(1/2);
f(2,1)=1/2;

t=zeros(vdim2,vdim2);
t(1,1)=Kinv(1,1)*(x1 - x2)^2 + (2*Kinv(1,2)*(x1 - x2) + Kinv(2,2)*(y1 - y2))*(y1 - y2);
t(2,1)=(Kinv(1,1)*(x1 - x2)*(x2 - x3)) + Kinv(2,2)*(y1 - y2)*(y2 - y3) + ...
         Kinv(1,2)*(x3*(-y1 + y2) + x1*(y2 - y3) + x2*(y1 - 2*y2 + y3));
t(1,2)=t(2,1);
t(2,2)=Kinv(1,1)*(x2 - x3)^2 + (2*Kinv(1,2)*(x2 - x3) + Kinv(2,2)*(y2 - y3))*(y2 - y3);
t=t/SE_denom_zz;

% CC.D en el nodo S-E, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
    % frontera Este
        Pg1u_L7=int_dir_simpson(i,j,i,j+1,tt,1);
        Pg2u_L7=int_dir_simpson(i,j,i,j+1,tt,2);
    % frontera Sur
        Pg1u_L3=int_dir_simpson(i-1,j,i,j,tt,1);
        Pg2u_L3=int_dir_simpson(i-1,j,i,j,tt,2);
    Pgu=(1/2)*[Pg1u_L7;Pg2u_L7;-Pg1u_L3;-Pg2u_L3];
    
localv_u=(b'-((b'*(a\c))*((c'*(a\c))\c')))*(a\Pgu);
gDu(ind1u:ind2u)=gDu(ind1u:ind2u)+localv_u;

% CC.D en el nodo S-E, fórmula de cuadratura para el término Pg de p en:
    % frontera Este
        Pgp_L7=int_dir_simpson(i,j,i,j+1,tt,3);
    % frontera Sur
        Pgp_L3=int_dir_simpson(i-1,j,i,j,tt,3);
    Pgp=(1/2)*[-Pgp_L7;Pgp_L3];
    
localv_p=((delta_t*f')*(t\Pgp))+((-d'+((d'*(a\c))*((c'*(a\c))\c')))*(a\Pgu));
gDp(ind1p)=gDp(ind1p)+localv_p;

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
    vdim2=vdim/2;
    
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

    JE1r4=abs(x3*y1 - x4*y1 - x1*y3 + x4*y3 + x1*y4 - x3*y4);
    W_denom_E1_ss=coef_denom_elas*JE1r4;
    W_denom_E1_zz=4*JE1r4;
    JE2r1=abs(x6*(-y3 + y4) + x4*(y3 - y6) + x3*(-y4 + y6));
    W_denom_E2_ss=coef_denom_elas*JE2r1;
    W_denom_E2_zz=4*JE2r1;
    
    a=zeros(vdim,vdim);    
    a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/W_denom_E1_ss;
    a(2,1)=-(lambda*(x3 - x4)*(y3 - y4))/W_denom_E1_ss;
    a(3,1)=(-2*mu*((x1 - x4)*(x3 - x4) + (y1 - y4)*(y3 - y4)) + lambda*((x1 - x4)*(-x3 + x4) + 2*(y1 - y4)*(-y3 + y4)))/W_denom_E1_ss;
    a(4,1)=(lambda*(x3 - x4)*(y1 - y4))/W_denom_E1_ss;
%     a(5,1)=0;
%     a(6,1)=0;
    a(1,2)=a(2,1);
    a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/W_denom_E1_ss;
    a(3,2)=(lambda*(x1 - x4)*(y3 - y4))/W_denom_E1_ss;
    a(4,2)=(-2*mu*((x1 - x4)*(x3 - x4) + (y1 - y4)*(y3 - y4)) + lambda*(2*(x1 - x4)*(-x3 + x4) + (y1 - y4)*(-y3 + y4)))/W_denom_E1_ss;
%     a(5,2)=0;
%     a(6,2)=0;
    a(1,3)=a(3,1);
    a(2,3)=a(3,2);
    a(3,3)=((lambda + 2*mu)*(x1 - x4)^2 + 2*(lambda + mu)*(y1 - y4)^2)/W_denom_E1_ss +...
           ((lambda + 2*mu)*(x4 - x6)^2 + 2*(lambda + mu)*(y4 - y6)^2)/W_denom_E2_ss;
    a(4,3)=-(lambda*(x1 - x4)*(y1 - y4))/W_denom_E1_ss - ...
           (lambda*(x4 - x6)*(y4 - y6))/W_denom_E2_ss;
    a(5,3)=-((lambda + 2*mu)*(x3 - x4)*(x4 - x6) + 2*(lambda + mu)*(y3 - y4)*(y4 - y6))/W_denom_E2_ss;
    a(6,3)=(lambda*(x4 - x6)*(y3 - y4))/W_denom_E2_ss;
    a(1,4)= a(4,1);
    a(2,4)= a(4,2);
    a(3,4)= a(4,3);
    a(4,4)=(2*(lambda + mu)*(x1 - x4)^2 + (lambda + 2*mu)*(y1 - y4)^2)/W_denom_E1_ss + ...
           (2*(lambda + mu)*(x4 - x6)^2 + (lambda + 2*mu)*(y4 - y6)^2)/W_denom_E2_ss;
    a(5,4)=(lambda*(x3 - x4)*(y4 - y6))/W_denom_E2_ss;
    a(6,4)=-(2*(lambda + mu)*(x3 - x4)*(x4 - x6) + (lambda + 2*mu)*(y3 - y4)*(y4 - y6))/W_denom_E2_ss;
%     a(1,5)=0;
%     a(2,5)=0;
    a(3,5)=a(5,3);
    a(4,5)=a(5,4);
    a(5,5)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/W_denom_E2_ss;
    a(6,5)=-(lambda*(x3 - x4)*(y3 - y4))/W_denom_E2_ss;
%     a(1,6)=0;
%     a(2,6)=0;
    a(3,6)=a(6,3);
    a(4,6)=a(6,4);
    a(5,6)=a(6,5);
    a(6,6)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/W_denom_E2_ss;
    
    b=zeros(vdim,4);
    b(1,1)=-(1/2);
%     b(1,2)=0;
%     b(2,1)=0;
    b(2,2)=-(1/2);
    b(3,1)=1/2;
%     b(3,2)=0;
    b(3,3)=-(1/2);
%     b(3,4)=0;
%     b(4,1)=0;
    b(4,2)=1/2;
%     b(4,3)=0;
    b(4,4)=-(1/2);
    b(5,3)=-(1/2);
%     b(5,4)=0;
%     b(6,3)=0;
    b(6,4)=-(1/2);
    
    c=(1/4)*[(y3 - y4),(-x3 + x4),(-y1 + y6),(x1 - x6),(y3 - y4),(-x3 + x4)]';
    
    d=zeros(vdim,2);
    d(1,1)=(x3 - x4);
    d(2,1)=(y3 - y4);
    d(3,1)=(-x1 + x4);
    d(4,1)=(-y1 + y4);
%     d(5,1)=0;
%     d(6,1)=0;
%     d(1,2)=0;
%     d(2,2)=0;
    d(3,2)=(-x4 + x6);
    d(4,2)=(-y4 + y6);
    d(5,2)=(x3 - x4);
    d(6,2)=(y3 - y4);
    d=coef_d*d;
    
    f=zeros(vdim2,2);
    f(1,1)=1/2;
    f(2,1)=-(1/2);
%     f(3,1)=0;
%     f(1,2)=0;
    f(2,2)=1/2;
    f(3,2)=1/2;
    
    t=zeros(vdim2,vdim2);  
    t(1,1)=(Kinv(1,1)*(x3 - x4)^2 + (2*Kinv(1,2)*(x3 - x4) + Kinv(2,2)*(y3 - y4))*(y3 - y4))/W_denom_E1_zz;
    t(2,1)=(Kinv(1,1)*(x1 - x4)*(-x3 + x4) + Kinv(2,2)*(y1 - y4)*(-y3 + y4) + ...
            Kinv(1,2)*(x4*(y1 + y3 - 2*y4) + x3*(-y1 + y4) + x1*(-y3 + y4)))/W_denom_E1_zz;
%     t(3,1)=0;
    t(1,2)=t(2,1);
    t(2,2)=(Kinv(1,1)*(x1 - x4)^2 + (2*Kinv(1,2)*(x1 - x4) + Kinv(2,2)*(y1 - y4))*(y1 - y4))/W_denom_E1_zz +...
           (Kinv(1,1)*(x4 - x6)^2 + (2*Kinv(1,2)*(x4 - x6) + Kinv(2,2)*(y4 - y6))*(y4 - y6))/W_denom_E2_zz;
    t(3,2)=((x4 - x6)*(Kinv(1,1)*(-x3 + x4) + Kinv(1,2)*(-y3 + y4)) + ...
            (Kinv(1,2)*(x3 - x4) + Kinv(2,2)*(y3 - y4))*(-y4 + y6))/W_denom_E2_zz;
%     t(1,3)=0;
    t(2,3)=t(3,2);
    t(3,3)=(Kinv(1,1)*(x3 - x4)^2 + (2*Kinv(1,2)*(x3 - x4) + Kinv(2,2)*(y3 - y4))*(y3 - y4))/W_denom_E2_zz;
    
    % CC.D en el nodo Oeste, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
        % frontera Oeste 1
            Pg1u_L4=int_dir_simpson(i,j,i,j-1,tt,1);
            Pg2u_L4=int_dir_simpson(i,j,i,j-1,tt,2);
        % frontera Oeste 2
            Pg1u_L11=int_dir_simpson(i,j,i,j+1,tt,1);
            Pg2u_L11=int_dir_simpson(i,j,i,j+1,tt,2);
        Pgu=(1/2)*[-Pg1u_L4;-Pg2u_L4;0;0;-Pg1u_L11;-Pg2u_L11];
        
    localv_u=(b'-((b'*(a\c))*((c'*(a\c))\c')))*(a\Pgu);
    gDu(ind1u:ind2u)=gDu(ind1u:ind2u)+localv_u(1:2);
    gDu(ind3u:ind4u)=gDu(ind3u:ind4u)+localv_u(3:4);
    
    % CC.D en el nodo Oeste, fórmula de cuadratura para el término Pg de p en:
        % frontera Oeste 1
            Pgp_L4=int_dir_simpson(i,j,i,j-1,tt,3);
        % frontera Oeste 2
            Pgp_L11=int_dir_simpson(i,j,i,j+1,tt,3);
        Pgp=(1/2)*[Pgp_L4;0;Pgp_L11];
        
    localv_p=((delta_t*f')*(t\Pgp))+((-d'+((d'*(a\c))*((c'*(a\c))\c')))*(a\Pgu));
    gDp(ind1p)=gDp(ind1p)+localv_p(1);
    gDp(ind2p)=gDp(ind2p)+localv_p(2);
    
    % Central nodes: No hay CCD
    
    % East nodes
    i=N+1;
    ind2u=(i-1+(j-2)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i-1+(j-1)*N)*2;
    ind3u=ind4u-1;
    ind1p=ind2u/2;
    ind2p=ind4u/2;
    
%     vdim=6;
%     vdim2=vdim/2;

    x1=x(i-1,j-1);
    y1=y(i-1,j-1);
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
    
    JE1r3=abs(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3) - ...
          abs(x2*y1 - x4*y1 - x1*y2 + x4*y2 + x1*y4 - x2*y4) +...
          abs(x3*y1 - x4*y1 - x1*y3 + x4*y3 + x1*y4 - x3*y4);
    E_denom_E1_ss=coef_denom_elas*JE1r3;
    E_denom_E1_zz=4*JE1r3;
    JE2r2=abs(x5*(-y3 + y4) + x4*(y3 - y5) + x3*(-y4 + y5));
    E_denom_E2_ss=coef_denom_elas*JE2r2;
    E_denom_E2_zz=4*JE2r2;
    
    a=zeros(vdim,vdim);
    a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/E_denom_E1_ss;
    a(2,1)=-(lambda*(x3 - x4)*(y3 - y4))/E_denom_E1_ss;
%     a(3,1)=0;
%     a(4,1)=0;
    a(5,1)=-((lambda + 2*mu)*(x2 - x3)*(x3 - x4) + 2*(lambda + mu)*(y2 - y3)*(y3 - y4))/E_denom_E1_ss;
    a(6,1)=(lambda*(x3 - x4)*(y2 - y3))/E_denom_E1_ss;
    a(1,2)=a(2,1);
    a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/E_denom_E1_ss;
%     a(3,2)=0;
%     a(4,2)=0;
    a(5,2)=(lambda*(x2 - x3)*(y3 - y4))/E_denom_E1_ss;
    a(6,2)=-(2*(lambda + mu)*(x2 - x3)*(x3 - x4) + (lambda + 2*mu)*(y2 - y3)*(y3 - y4))/E_denom_E1_ss;
%     a(1,3)=0;
%     a(2,3)=0;
    a(3,3)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/E_denom_E2_ss;
    a(4,3)=-(lambda*(x3 - x4)*(y3 - y4))/E_denom_E2_ss;
    a(5,3)=-((lambda + 2*mu)*(x3 - x4)*(x3 - x5) + 2*(lambda + mu)*(y3 - y4)*(y3 - y5))/E_denom_E2_ss;
    a(6,3)=(lambda*(x3 - x4)*(y3 - y5))/E_denom_E2_ss;
%     a(1,4)=0;
%     a(2,4)=0;
    a(3,4)=a(4,3);
    a(4,4)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/E_denom_E2_ss;
    a(5,4)=(lambda*(x3 - x5)*(y3 - y4))/E_denom_E2_ss;
    a(6,4)=-(2*(lambda + mu)*(x3 - x4)*(x3 - x5) + (lambda + 2*mu)*(y3 - y4)*(y3 - y5))/E_denom_E2_ss;
    a(1,5)=a(5,1);
    a(2,5)=a(5,2);
    a(3,5)=a(5,3);
    a(4,5)=a(5,4);
    a(5,5)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/E_denom_E1_ss +...
           ((lambda + 2*mu)*(x3 - x5)^2 + 2*(lambda + mu)*(y3 - y5)^2)/E_denom_E2_ss;
    a(6,5)=-(lambda*(x2 - x3)*(y2 - y3))/E_denom_E1_ss -...
           (lambda*(x3 - x5)*(y3 - y5))/E_denom_E2_ss;
    a(1,6)=a(6,1);
    a(2,6)=a(6,2);
    a(3,6)=a(6,3);
    a(4,6)=a(6,4);
    a(5,6)=a(6,5);
    a(6,6)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/E_denom_E1_ss +...
           (2*(lambda + mu)*(x3 - x5)^2 + (lambda + 2*mu)*(y3 - y5)^2)/E_denom_E2_ss;
    
    b=zeros(vdim,4);
    b(1,1)=1/2;
%     b(1,2)=0;
%     b(2,1)=0;
    b(2,2)=1/2;
    b(3,3)=1/2;
%     b(3,4)=0;
%     b(4,3)=0;
    b(4,4)=1/2;
    b(5,1)=1/2;
%     b(5,2)=0;
    b(5,3)=-(1/2);
%     b(5,4)=0;
%     b(6,1)=0;
    b(6,2)=1/2;
%     b(6,3)=0;
    b(6,4)=-(1/2);
    
    c=(1/4)*[(y3 - y4),(-x3 + x4),(y3 - y4),(-x3 + x4),(-y2 + y5),(x2 - x5)]';
    
    d=zeros(vdim,2);
    d(1,1)=(x3 - x4);
    d(2,1)=(y3 - y4);
%     d(3,1)=0;
%     d(4,1)=0;
    d(5,1)=(-x2 + x3);
    d(6,1)=(-y2 + y3);
%     d(1,2)=0;
%     d(2,2)=0;
    d(3,2)=(x3 - x4);
    d(4,2)=(y3 - y4);
    d(5,2)=(-x3 + x5);
    d(6,2)=(-y3 + y5);
    d=coef_d*d;
    
    f=zeros(vdim2,2);
    f(1,1)=-(1/2);
%     f(2,1)=0;
    f(3,1)=-(1/2);
%     f(1,2)=0;
    f(2,2)=-(1/2);
    f(3,2)=1/2;
    
    t=zeros(vdim2,vdim2);
    t(1,1)=(Kinv(1,1)*(x3 - x4)^2 + (2*Kinv(1,2)*(x3 - x4) + Kinv(2,2)*(y3 - y4))*(y3 - y4))/E_denom_E1_zz;
%     t(2,1)=0;
    t(3,1)=-(Kinv(1,1)*(x2 - x3)*(x3 - x4) + Kinv(2,2)*(y2 - y3)*(y3 - y4) + ...
             Kinv(1,2)*(x4*(-y2 + y3) + x2*(y3 - y4) + x3*(y2 - 2*y3 + y4)))/E_denom_E1_zz;
%     t(1,2)=0;
    t(2,2)=(Kinv(1,1)*(x3 - x4)^2 + (2*Kinv(1,2)*(x3 - x4) + Kinv(2,2)*(y3 - y4))*(y3 - y4))/E_denom_E2_zz;
    t(3,2)=-(Kinv(1,1)*(x3 - x4)*(x3 - x5) + Kinv(2,2)*(y3 - y4)*(y3 - y5) - ...
             Kinv(1,2)*(x5*(y3 - y4) + x4*(y3 - y5) + x3*(-2*y3 + y4 + y5)))/E_denom_E2_zz;
    t(1,3)=t(3,1);
    t(2,3)=t(3,2);
    t(3,3)=(Kinv(1,1)*(x2 - x3)^2 + (2*Kinv(1,2)*(x2 - x3) + Kinv(2,2)*(y2 - y3))*(y2 - y3))/E_denom_E1_zz +...
           (Kinv(1,1)*(x3 - x5)^2 + (2*Kinv(1,2)*(x3 - x5) + Kinv(2,2)*(y3 - y5))*(y3 - y5))/E_denom_E2_zz;
    
    % CC.D en el nodo Este, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
        % frontera Este 1
            Pg1u_L7=int_dir_simpson(i,j,i,j-1,tt,1);
            Pg2u_L7=int_dir_simpson(i,j,i,j-1,tt,2);
        % frontera Este 2
            Pg1u_L14=int_dir_simpson(i,j,i,j+1,tt,1);
            Pg2u_L14=int_dir_simpson(i,j,i,j+1,tt,2);
        Pgu=(1/2)*[Pg1u_L7;Pg2u_L7;Pg1u_L14;Pg2u_L14;0;0];
    
       localv_u=(b'-((b'*(a\c))*((c'*(a\c))\c')))*(a\Pgu);
       gDu(ind1u:ind2u)=gDu(ind1u:ind2u)+localv_u(1:2);
       gDu(ind3u:ind4u)=gDu(ind3u:ind4u)+localv_u(3:4);
       
       % CC.D en el nodo Este, fórmula de cuadratura para el término Pg de p en:
        % frontera Este 1
            Pgp_L7=int_dir_simpson(i,j,i,j-1,tt,3);
        % frontera Este 2
            Pgp_L14=int_dir_simpson(i,j,i,j+1,tt,3);
        Pgp=-(1/2)*[Pgp_L7;Pgp_L14;0];
    
       localv_p=((delta_t*f')*(t\Pgp))+((-d'+((d'*(a\c))*((c'*(a\c))\c')))*(a\Pgu));
       gDp(ind1p)=gDp(ind1p)+localv_p(1);
       gDp(ind2p)=gDp(ind2p)+localv_p(2);
end

% North-West corner node
i=1;
j=N+1;
ind2u=(i+(j-2)*N)*2;
ind1u=ind2u-1;
ind1p=ind2u/2;
vdim=4;
vdim2=vdim/2;

x1=x(i,j-1);
y1=y(i,j-1);           
% x2=x(i+1,j-1);
% y2=y(i+1,j-1);
x3=x(i+1,j);
y3=y(i+1,j);
x4=x(i,j);
y4=y(i,j);

Jer4=abs(x3*y1 - x4*y1 - x1*y3 + x4*y3 + x1*y4 - x3*y4);
NW_denom_ss=coef_denom_elas*Jer4;
NW_denom_zz=4*Jer4;

a=zeros(vdim,vdim);
a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2);
a(2,1)=-(lambda*(x3 - x4)*(y3 - y4));
a(3,1)=(-2*mu*((x1 - x4)*(x3 - x4) + (y1 - y4)*(y3 - y4)) + lambda*((x1 - x4)*(-x3 + x4) + 2*(y1 - y4)*(-y3 + y4)));
a(4,1)=(lambda*(x3 - x4)*(y1 - y4));
a(1,2)=a(2,1);
a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2);
a(3,2)=(lambda*(x1 - x4)*(y3 - y4));
a(4,2)=(-2*mu*((x1 - x4)*(x3 - x4) + (y1 - y4)*(y3 - y4)) + lambda*(2*(x1 - x4)*(-x3 + x4) + (y1 - y4)*(-y3 + y4)));
a(1,3)=a(3,1);
a(2,3)=a(3,2);
a(3,3)=((lambda + 2*mu)*(x1 - x4)^2 + 2*(lambda + mu)*(y1 - y4)^2);
a(4,3)=-(lambda*(x1 - x4)*(y1 - y4));
a(1,4)=a(4,1);
a(2,4)=a(4,2);
a(3,4)=a(4,3);
a(4,4)=(2*(lambda + mu)*(x1 - x4)^2 + (lambda + 2*mu)*(y1 - y4)^2);
a=a/NW_denom_ss;

b=zeros(vdim,2);
b(1,1)=-(1/2);
% b(1,2)=20;
% b(2,1)=0;
b(2,2)=-(1/2);
b(3,1)=1/2;
% b(3,2)=0;
% b(4,1)=0;
b(4,2)=1/2;

c=(1/4)*[(y3 - y4),(-x3 + x4),(-y1 + y4),(x1 - x4)]';

d=zeros(vdim,1);
d(1,1)=(x3 - x4);
d(2,1)=(y3 - y4);
d(3,1)=(-x1 + x4);
d(4,1)=(-y1 + y4);
d=coef_d*d;

f=zeros(vdim2,1);
f(1,1)=1/2;
f(2,1)=-(1/2);

t=zeros(vdim2,vdim2);
t(1,1)=Kinv(1,1)*(x3 - x4)^2 + (2*Kinv(1,2)*(x3 - x4) + Kinv(2,2)*(y3 - y4))*(y3 - y4);
t(2,1)=Kinv(1,1)*(x1 - x4)*(-x3 + x4) + Kinv(2,2)*(y1 - y4)*(-y3 + y4) + ...
       Kinv(1,2)*(x4*(y1 + y3 - 2*y4) + x3*(-y1 + y4) + x1*(-y3 + y4));
t(1,2)=t(2,1);
t(2,2)=Kinv(1,1)*(x1 - x4)^2 + (2*Kinv(1,2)*(x1 - x4) + Kinv(2,2)*(y1 - y4))*(y1 - y4);
t=t/NW_denom_zz;

% CC.D en el nodo N-W, fórmula de cuadratura para los término Pg1 y Pg2 en:
    % frontera Oeste
        Pg1u_L18=int_dir_simpson(i,j,i,j-1,tt,1);
        Pg2u_L18=int_dir_simpson(i,j,i,j-1,tt,2);
    % frontera Norte
        Pg1u_L22=int_dir_simpson(i,j,i+1,j,tt,1);
        Pg2u_L22=int_dir_simpson(i,j,i+1,j,tt,2);
    Pgu=(1/2)*[-Pg1u_L18;-Pg2u_L18;Pg1u_L22;Pg2u_L22];
    
localv_u=(b'-((b'*(a\c))*((c'*(a\c))\c')))*(a\Pgu);
gDu(ind1u:ind2u)=gDu(ind1u:ind2u)+localv_u;

% CC.D en el nodo N-W, fórmula de cuadratura para el término Pg de p en:
    % frontera Oeste
        Pgp_L18=int_dir_simpson(i,j,i,j-1,tt,3);
    % frontera Norte
        Pgp_L22=int_dir_simpson(i,j,i+1,j,tt,3);
    Pgp=(1/2)*[Pgp_L18;-Pgp_L22];
    
localv_p=((delta_t*f')*(t\Pgp))+((-d'+((d'*(a\c))*((c'*(a\c))\c')))*(a\Pgu));
gDp(ind1p)=gDp(ind1p)+localv_p;

% North nodes (j=N+1)
vdim=6;
vdim2=vdim/2;

% Matrix b and f are the same for every North node
b=zeros(vdim,4);
b(1,1)=1/2;
% b(1,2)=0;
b(1,3)=-(1/2);
% b(1,4)=0;
% b(2,1)=0;
b(2,2)=1/2;
% b(2,3)=0;
b(2,4)=-(1/2);
b(3,3)=1/2;
% b(3,4)=0;
% b(4,3)=0;
b(4,4)=1/2;
b(5,1)=1/2;
% b(5,2)=0;
% b(6,1)=0;
b(6,2)=1/2;

f=zeros(vdim2,2);
f(1,1)=-(1/2);
% f(2,1)=0;
f(3,1)=-(1/2);
f(1,2)=1/2;
f(2,2)=-(1/2);
% f(3,2)=0;

a=zeros(vdim,vdim);
d=zeros(vdim,2);
t=zeros(vdim2,vdim2);

for i=2:N
    ind2u=(i-1+(j-2)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i+(j-2)*N)*2;
%     ind3u=ind4u-1;
    ind1p=ind2u/2;
    ind2p=ind4u/2;
    
    x1=x(i-1,j-1);
    y1=y(i-1,j-1);
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
    
    JE1r3=abs(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3) - ...
          abs(x2*y1 - x4*y1 - x1*y2 + x4*y2 + x1*y4 - x2*y4) + ...
          abs(x3*y1 - x4*y1 - x1*y3 + x4*y3 + x1*y4 - x3*y4);
    N_denom_E1_ss=coef_denom_elas*JE1r3;
    N_denom_E1_zz=4*JE1r3;
    JE2r4=abs(x6*(-y2 + y3) + x3*(y2 - y6) + x2*(-y3 + y6));
    N_denom_E2_ss=coef_denom_elas*JE2r4;
    N_denom_E2_zz=4*JE2r4;
    
    a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/N_denom_E1_ss +...
           ((lambda + 2*mu)*(x3 - x6)^2 + 2*(lambda + mu)*(y3 - y6)^2)/N_denom_E2_ss;
    a(2,1)=-(lambda*(x3 - x4)*(y3 - y4))/N_denom_E1_ss -...
           (lambda*(x3 - x6)*(y3 - y6))/N_denom_E2_ss;
    a(3,1)=((lambda + 2*mu)*(x2 - x3)*(x3 - x6) + 2*(lambda + mu)*(y2 - y3)*(y3 - y6))/N_denom_E2_ss;
    a(4,1)=-(lambda*(x3 - x6)*(y2 - y3))/N_denom_E2_ss;
    a(5,1)=-((lambda + 2*mu)*(x2 - x3)*(x3 - x4) + 2*(lambda + mu)*(y2 - y3)*(y3 - y4))/N_denom_E1_ss;
    a(6,1)=(lambda*(x3 - x4)*(y2 - y3))/N_denom_E1_ss;
    a(1,2)=a(2,1);
    a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/N_denom_E1_ss +...
           (2*(lambda + mu)*(x3 - x6)^2 + (lambda + 2*mu)*(y3 - y6)^2)/N_denom_E2_ss;
    a(3,2)=-(lambda*(x2 - x3)*(y3 - y6))/N_denom_E2_ss;
    a(4,2)=(2*(lambda + mu)*(x2 - x3)*(x3 - x6) + (lambda + 2*mu)*(y2 - y3)*(y3 - y6))/N_denom_E2_ss;
    a(5,2)=(lambda*(x2 - x3)*(y3 - y4))/N_denom_E1_ss;
    a(6,2)=-(2*(lambda + mu)*(x2 - x3)*(x3 - x4) + (lambda + 2*mu)*(y2 - y3)*(y3 - y4))/N_denom_E1_ss;
    a(1,3)=a(3,1);
    a(2,3)=a(3,2);
    a(3,3)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/N_denom_E2_ss;
    a(4,3)=-(lambda*(x2 - x3)*(y2 - y3))/N_denom_E2_ss;
%     a(5,3)=0;
%     a(6,3)=0;
    a(1,4)=a(4,1);
    a(2,4)=a(4,2);
    a(3,4)=a(4,3);
    a(4,4)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/N_denom_E2_ss;
%     a(5,4)=0;
%     a(6,4)=0;
    a(1,5)=a(5,1);
    a(2,5)=a(5,2);
%     a(3,5)=0;
%     a(4,5)=0;
    a(5,5)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/N_denom_E1_ss;
    a(6,5)=-(lambda*(x2 - x3)*(y2 - y3))/N_denom_E1_ss;
    a(1,6)=a(6,1);
    a(2,6)=a(6,2);
%     a(3,6)=0;
%     a(4,6)=0;
    a(5,6)=a(6,5);
    a(6,6)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/N_denom_E1_ss;
    
    c=(1/4)*[(-y4 + y6),(x4 - x6),(-y2 + y3),(x2 - x3),(-y2 + y3),(x2 - x3)]';
    
    d(1,1)=(x3 - x4);
    d(2,1)=(y3 - y4);
%     d(3,1)=0;
%     d(4,1)=0;
    d(5,1)=(-x2 + x3);
    d(6,1)=(-y2 + y3);
    d(1,2)=(-x3 + x6);
    d(2,2)=(-y3 + y6);
    d(3,2)=(-x2 + x3);
    d(4,2)=(-y2 + y3);
%     d(5,2)=0;
%     d(6,2)=0;
    d=coef_d*d;

    t(1,1)=(Kinv(1,1)*(x3 - x4)^2 + (2*Kinv(1,2)*(x3 - x4) + Kinv(2,2)*(y3 - y4))*(y3 - y4))/N_denom_E1_zz +...
           (Kinv(1,1)*(x3 - x6)^2 + (2*Kinv(1,2)*(x3 - x6) + Kinv(2,2)*(y3 - y6))*(y3 - y6))/N_denom_E2_zz;
    t(2,1)=(Kinv(1,1)*(x2 - x3)*(x3 - x6) + Kinv(2,2)*(y2 - y3)*(y3 - y6) + ...
            Kinv(1,2)*(x6*(-y2 + y3) + x2*(y3 - y6) + x3*(y2 - 2*y3 + y6)))/N_denom_E2_zz;
    t(3,1)=-(Kinv(1,1)*(x2 - x3)*(x3 - x4) + Kinv(2,2)*(y2 - y3)*(y3 - y4) + ...
             Kinv(1,2)*(x4*(-y2 + y3) + x2*(y3 - y4) + x3*(y2 - 2*y3 + y4)))/N_denom_E1_zz;
    t(1,2)=t(2,1);
    t(2,2)=(Kinv(1,1)*(x2 - x3)^2 + (2*Kinv(1,2)*(x2 - x3) + Kinv(2,2)*(y2 - y3))*(y2 - y3))/N_denom_E2_zz;
%     t(3,2)=0;
    t(1,3)=t(3,1);
%     t(2,3)=0;
    t(3,3)=(Kinv(1,1)*(x2 - x3)^2 + (2*Kinv(1,2)*(x2 - x3) + Kinv(2,2)*(y2 - y3))*(y2 - y3))/N_denom_E1_zz;
    
    % CC.D en el nodo Norte, fórmula de cuadratura para los términos Pg1 y Pg2 de u en:
        % frontera Norte 1
            Pg1u_L23=int_dir_simpson(i,j,i+1,j,tt,1);
            Pg2u_L23=int_dir_simpson(i,j,i+1,j,tt,2);
        % frontera Norte 2
            Pg1u_L22=int_dir_simpson(i-1,j,i,j,tt,1);
            Pg2u_L22=int_dir_simpson(i-1,j,i,j,tt,2);
        Pgu=(1/2)*[0;0;Pg1u_L23;Pg2u_L23;Pg1u_L22;Pg2u_L22];
    
    localv_u=(b'-((b'*(a\c))*((c'*(a\c))\c')))*(a\Pgu);
    gDu(ind1u:ind4u)=gDu(ind1u:ind4u)+localv_u;
    
    % CC.D en el nodo Norte, fórmula de cuadratura para el término Pg de p en:
        % frontera Norte 1
            Pgp_L23=int_dir_simpson(i,j,i+1,j,tt,3);
        % frontera Norte 2
            Pgp_L22=int_dir_simpson(i-1,j,i,j,tt,3);
        Pgp=-(1/2)*[0;Pgp_L23;Pgp_L22];
    
    localv_p=((delta_t*f')*(t\Pgp))+((-d'+((d'*(a\c))*((c'*(a\c))\c')))*(a\Pgu));
    gDp(ind1p:ind2p)=gDp(ind1p:ind2p)+localv_p;
end

% North-East corner node (j=N+1)
i=N+1;
ind2u=(i-1+(j-2)*N)*2;
ind1u=ind2u-1;
ind1p=ind2u/2;
vdim=4;
vdim2=vdim/2;
          
x1=x(i-1,j-1);  
y1=y(i-1,j-1);            
x2=x(i,j-1);
y2=y(i,j-1);
x3=x(i,j);
y3=y(i,j);
x4=x(i-1,j);
y4=y(i-1,j);

Jer3=abs(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3) -...
    abs(x2*y1 - x4*y1 - x1*y2 + x4*y2 + x1*y4 - x2*y4) +...
    abs(x3*y1 - x4*y1 - x1*y3 + x4*y3 + x1*y4 - x3*y4);
NE_denom_ss=coef_denom_elas*Jer3;
NE_denom_zz=4*Jer3;

a=zeros(vdim,vdim);
a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2);
a(2,1)=-(lambda*(x3 - x4)*(y3 - y4));
a(3,1)=-((lambda + 2*mu)*(x2 - x3)*(x3 - x4) + 2*(lambda + mu)*(y2 - y3)*(y3 - y4));
a(4,1)=(lambda*(x3 - x4)*(y2 - y3));
a(1,2)=a(2,1);
a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2);
a(3,2)=(lambda*(x2 - x3)*(y3 - y4));
a(4,2)=-(2*(lambda + mu)*(x2 - x3)*(x3 - x4) + (lambda + 2*mu)*(y2 - y3)*(y3 - y4));
a(1,3)=a(3,1);
a(2,3)=a(3,2);
a(3,3)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2);
a(4,3)=-(lambda*(x2 - x3)*(y2 - y3));
a(1,4)=a(4,1);
a(2,4)=a(4,2);
a(3,4)=a(4,3);
a(4,4)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2);
a=a/NE_denom_ss;

b=zeros(vdim,2);
b(1,1)=1/2;
% b(1,2)=0;
% b(2,1)=0;
b(2,2)=1/2;
b(3,1)=1/2;
% b(3,2)=0;
% b(4,1)=0;
b(4,2)=1/2;

c=(1/4)*[(y3 - y4),(-x3 + x4),(-y2 + y3),(x2 - x3)]';

d=zeros(vdim,1);
d(1,1)=(x3 - x4);
d(2,1)=(y3 - y4);
d(3,1)=(-x2 + x3);
d(4,1)=(-y2 + y3);
d=coef_d*d;

f=zeros(vdim2,1);
f(1,1)=-(1/2);
f(2,1)=-(1/2);

t=zeros(vdim2,vdim2);
t(1,1)=Kinv(1,1)*(x3 - x4)^2 + (2*Kinv(1,2)*(x3 - x4) + Kinv(2,2)*(y3 - y4))*(y3 - y4);
t(2,1)=-(Kinv(1,1)*(x2 - x3)*(x3 - x4)) + Kinv(2,2)*(y2 - y3)*(y3 - y4) + ...
         Kinv(1,2)*(x4*(-y2 + y3) + x2*(y3 - y4) + x3*(y2 - 2*y3 + y4));
t(1,2)=t(2,1);
t(2,2)=Kinv(1,1)*(x2 - x3)^2 + (2*Kinv(1,2)*(x2 - x3) + Kinv(2,2)*(y2 - y3))*(y2 - y3);
t=t/NE_denom_zz;

% CC.D en el nodo N-E, fórmula de cuadratura para los términos Pg1 y Pg2 de u en:
    % frontera Este
        Pg1u_L21=int_dir_simpson(i,j,i,j-1,tt,1);
        Pg2u_L21=int_dir_simpson(i,j,i,j-1,tt,2);
    % frontera Norte
        Pg1u_L24=int_dir_simpson(i-1,j,i,j,tt,1);
        Pg2u_L24=int_dir_simpson(i-1,j,i,j,tt,2);
    Pgu=(1/2)*[Pg1u_L21;Pg2u_L21;Pg1u_L24;Pg2u_L24];
    
localv_u=(b'-((b'*(a\c))*((c'*(a\c))\c')))*(a\Pgu);
gDu(ind1u:ind2u)=gDu(ind1u:ind2u)+localv_u;

% CC.D en el nodo N-E, fórmula de cuadratura para el término Pg de p en:
    % frontera Este
        Pgp_L21=int_dir_simpson(i,j,i,j-1,tt,3);
    % frontera Norte
        Pgp_L24=int_dir_simpson(i-1,j,i,j,tt,3);
    Pgp=-(1/2)*[Pgp_L21;Pgp_L24];
    
localv_p=((delta_t*f')*(t\Pgp))+((-d'+((d'*(a\c))*((c'*(a\c))\c')))*(a\Pgu));
gDp(ind1p)=gDp(ind1p)+localv_p;
return
end
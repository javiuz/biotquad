function [z,zx,zy]=compute_fluxes(paux,tt)

global NN x y Kinv

N=NN;

z=zeros(4*N*(N+1),1);

% p=zeros(nx,ny);
% 
% for j=1:ny
%     for i=1:nx
%         p(i,j)=paux(i+(j-1)*ny);
%     end
% end

p=reshape(paux,N,N);

zx=zeros(N+1,N+1,2);
zy=zx;

% South-West corner node
i=1;
j=1;
vdim=2;
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

Jer1=abs(x2*y1 - x4*y1 - x1*y2 + x4*y2 + x1*y4 - x2*y4);
SW_denom_zz=Jer1*4;

f=zeros(vdim,1);
f(1,1)=1/2;
f(2,1)=1/2;

t=zeros(vdim,vdim);
t(1,1)=Kinv(1,1)*(x1 - x4)^2 + (2*Kinv(1,2)*(x1 - x4) + Kinv(2,2)*(y1 - y4))*(y1 - y4);
t(2,1)=(x1 - x4)*(Kinv(1,1)*(x1 - x2) + Kinv(1,2)*(y1 - y2)) + ...
       (Kinv(1,2)*(x1 - x2) + Kinv(2,2)*(y1 - y2))*(y1 - y4);
t(1,2)=t(2,1);
t(2,2)=Kinv(1,1)*(x1 - x2)^2 + (2*Kinv(1,2)*(x1 - x2) + Kinv(2,2)*(y1 - y2))*(y1 - y2);
t=t/SW_denom_zz;

% CC.D en el nodo S-W, fórmula de cuadratura para el término Pg de p en:
    % frontera Sur
        Pgp_L1=int_dir_simpson(i,j,i+1,j,tt,3);
    % frontera Oeste
        Pgp_L4=int_dir_simpson(i,j,i,j+1,tt,3); 
        Pgp=(1/2)*[Pgp_L1;Pgp_L4];

pp(1,1)=p(i,j);
zz=t\(Pgp-f*pp);

z(ld:vdim)=zz;
ld=vdim;

zy(i,j,1)=zz(1,1);
zx(i,j,2)=zz(2,1);

% South nodes (j=1)
vdim=3;
pdim=2;
pp=zeros(pdim,1);

% Matrix f is the same for every South node
f=zeros(vdim,2);
% f(1,1)=0;
f(2,1)=-(1/2);
f(3,1)=1/2;
f(1,2)=1/2;
f(2,2)=1/2;
% f(3,2)=0;

t=zeros(vdim,vdim);

for i=2:N
    
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
    S_denom_E1_zz=4*JE1r2;
    JE2r1=abs(x5*(-y2 + y3) + x3*(y2 - y5) + x2*(-y3 + y5));
    S_denom_E2_zz=4*JE2r1;

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
    
    % CC.D en los nodos Sur, fórmula de cuadratura para el término Pg de p en:
        % frontera Sur 1
            Pgp_L2=int_dir_simpson(i,j,i+1,j,tt,3); 
        % frontera Sur 2
            Pgp_L1=int_dir_simpson(i-1,j,i,j,tt,3); 
        Pgp=(1/2)*[Pgp_L2;0;Pgp_L1];
    
    pp(1,1)=p(i-1,j);
    pp(2,1)=p(i,j);
    zz=t\(Pgp-f*pp);
    
    z(ld+1:ld+vdim)=zz;
    ld=ld+vdim;

    zy(i,j,1)=zz(1,1);
    zx(i,j,2)=zz(2,1);
    zy(i,j,2)=zz(3,1);
end 
    
% South-East corner node (j=1)
i=N+1;
vdim=2;
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

Jer2=abs(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3);
SE_denom_zz=Jer2*4;

f=zeros(vdim,1);
f(1,1)=-(1/2);
f(2,1)=1/2;

t=zeros(vdim,vdim);
t(1,1)=Kinv(1,1)*(x1 - x2)^2 + (2*Kinv(1,2)*(x1 - x2) + Kinv(2,2)*(y1 - y2))*(y1 - y2);
t(2,1)=(Kinv(1,1)*(x1 - x2)*(x2 - x3)) + Kinv(2,2)*(y1 - y2)*(y2 - y3) + ...
         Kinv(1,2)*(x3*(-y1 + y2) + x1*(y2 - y3) + x2*(y1 - 2*y2 + y3));
t(1,2)=t(2,1);
t(2,2)=Kinv(1,1)*(x2 - x3)^2 + (2*Kinv(1,2)*(x2 - x3) + Kinv(2,2)*(y2 - y3))*(y2 - y3);
t=t/SE_denom_zz;

% CC.D en el nodo S-E, fórmula de cuadratura para el término Pg de p en:
    % frontera Este
        Pgp_L7=int_dir_simpson(i,j,i,j+1,tt,3);
    % frontera Sur
        Pgp_L3=int_dir_simpson(i-1,j,i,j,tt,3);
    Pgp=(1/2)*[-Pgp_L7;Pgp_L3];

pp(1,1)=p(i-1,j);
zz=t\(Pgp-f*pp);

z(ld+1:ld+vdim)=zz;
ld=ld+vdim;

zx(i,j,2)=zz(1,1);
zy(i,j,2)=zz(2,1);

for j=2:N

    % West nodes
    i=1;
    vdim=3;
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
    
    JE1r4=abs(x3*y1 - x4*y1 - x1*y3 + x4*y3 + x1*y4 - x3*y4);
    W_denom_E1_zz=4*JE1r4;
    JE2r1=abs(x6*(-y3 + y4) + x4*(y3 - y6) + x3*(-y4 + y6));
    W_denom_E2_zz=4*JE2r1;
    
    f=zeros(vdim,2);
    f(1,1)=1/2;
    f(2,1)=-(1/2);
%     f(3,1)=0;
%     f(1,2)=0;
    f(2,2)=1/2;
    f(3,2)=1/2;
    
    t=zeros(vdim,vdim);
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
    
    % CC.D en el nodo Oeste, fórmula de cuadratura para el término Pg de p en:
        % frontera Oeste 1
            Pgp_L4=int_dir_simpson(i,j,i,j-1,tt,3);
        % frontera Oeste 2
            Pgp_L11=int_dir_simpson(i,j,i,j+1,tt,3);
        Pgp=(1/2)*[Pgp_L4;0;Pgp_L11];
        
    pp(1,1)=p(i,j-1);
    pp(2,1)=p(i,j);
    zz=t\(Pgp-f*pp);
    
    z(ld+1:ld+vdim)=zz;
    ld=ld+vdim;
    
    zx(i,j,1)=zz(1,1);
    zy(i,j,1)=zz(2,1);
    zx(i,j,2)=zz(3,1);

    % Central nodes 
    vdim=4;
    pdim=4;
    pp=zeros(pdim,1);
    
    t=zeros(vdim,vdim);

    % Matrix f is the same for every central node
    f=zeros(vdim,vdim);
    f(1,1)=-(1/2);
%     f(2,1)=0;
%     f(3,1)=0;
    f(4,1)=-(1/2);
    f(1,2)=1/2;
    f(2,2)=-(1/2);
%     f(3,2)=0;
%     f(4,2)=0;
%     f(1,3)=0;
%     f(2,3)=0;
    f(3,3)=-(1/2);
    f(4,3)=1/2;
%     f(1,4)=0;
    f(2,4)=1/2;
    f(3,4)=1/2;
%     f(4,4)=0;
    
    for i=2:N
        
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
%     x7=x(i+1,j+1);
%     y7=y(i+1,j+1);
    x8=x(i,j+1);
    y8=y(i,j+1);
%     x9=x(i-1,j+1);
%     y9=y(i-1,j+1);

    JE1r3=abs(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3) - ...
          abs(x2*y1 - x4*y1 - x1*y2 + x4*y2 + x1*y4 - x2*y4) + ...
          abs(x3*y1 - x4*y1 - x1*y3 + x4*y3 + x1*y4 - x3*y4);
    I_denom_E1_zz=4*JE1r3;
    JE2r4=abs(x6*(-y2 + y3) + x3*(y2 - y6) + x2*(-y3 + y6));4;
    I_denom_E2_zz=4*JE2r4;
    JE3r1=abs(x6*y3 - x8*y3 - x3*y6 + x8*y6 + x3*y8 - x6*y8);
    I_denom_E3_zz=4*JE3r1;
    JE4r2=abs(x8*(-y3 + y4) + x4*(y3 - y8) + x3*(-y4 + y8));
    I_denom_E4_zz=4*JE4r2;

    t(1,1)=(Kinv(1,1)*(x3 - x4)^2 + (2*Kinv(1,2)*(x3 - x4) + Kinv(2,2)*(y3 - y4))*(y3 - y4))/I_denom_E1_zz +...
           (Kinv(1,1)*(x3 - x6)^2 + (2*Kinv(1,2)*(x3 - x6) + Kinv(2,2)*(y3 - y6))*(y3 - y6))/I_denom_E2_zz;
    t(2,1)=(Kinv(1,1)*(x2 - x3)*(x3 - x6) + Kinv(2,2)*(y2 - y3)*(y3 - y6) + ...
            Kinv(1,2)*(x6*(-y2 + y3) + x2*(y3 - y6) + x3*(y2 - 2*y3 + y6)))/I_denom_E2_zz;
%     t(3,1)=0;
    t(4,1)=-(Kinv(1,1)*(x2 - x3)*(x3 - x4) + Kinv(2,2)*(y2 - y3)*(y3 - y4) + ...
             Kinv(1,2)*(x4*(-y2 + y3) + x2*(y3 - y4) + x3*(y2 - 2*y3 + y4)))/I_denom_E1_zz;
    t(1,2)=t(2,1);
    t(2,2)=(Kinv(1,1)*(x2 - x3)^2 + (2*Kinv(1,2)*(x2 - x3) + Kinv(2,2)*(y2 - y3))*(y2 - y3))/I_denom_E2_zz +...
           (Kinv(1,1)*(x3 - x8)^2 + (2*Kinv(1,2)*(x3 - x8) + Kinv(2,2)*(y3 - y8))*(y3 - y8))/I_denom_E3_zz;
    t(3,2)=((x3 - x8)*(Kinv(1,1)*(x3 - x6) + Kinv(1,2)*(y3 - y6)) + ...
           (Kinv(1,2)*(x3 - x6) + Kinv(2,2)*(y3 - y6))*(y3 - y8))/I_denom_E3_zz;
%     t(4,2)=0;
%     t(1,3)=0;
    t(2,3)=t(3,2);
    t(3,3)=(Kinv(1,1)*(x3 - x6)^2 + (2*Kinv(1,2)*(x3 - x6) + Kinv(2,2)*(y3 - y6))*(y3 - y6))/I_denom_E3_zz +...
           (Kinv(1,1)*(x3 - x4)^2 + (2*Kinv(1,2)*(x3 - x4) + Kinv(2,2)*(y3 - y4))*(y3 - y4))/I_denom_E4_zz;
    t(4,3)=-(Kinv(1,1)*(x3 - x4)*(x3 - x8) + Kinv(2,2)*(y3 - y4)*(y3 - y8) - ...
             Kinv(1,2)*(x8*(y3 - y4) + x4*(y3 - y8) + x3*(-2*y3 + y4 + y8)))/I_denom_E4_zz;
    t(1,4)=t(4,1);
%     t(2,4)=0;
    t(3,4)=t(4,3);
    t(4,4)=(Kinv(1,1)*(x2 - x3)^2 + (2*Kinv(1,2)*(x2 - x3) + Kinv(2,2)*(y2 - y3))*(y2 - y3))/I_denom_E1_zz +...
           (Kinv(1,1)*(x3 - x8)^2 + (2*Kinv(1,2)*(x3 - x8) + Kinv(2,2)*(y3 - y8))*(y3 - y8))/I_denom_E4_zz;
        
        pp(1,1)=p(i-1,j-1);
        pp(2,1)=p(i,j-1);
        pp(3,1)=p(i-1,j);
        pp(4,1)=p(i,j);
        zz=-t\(f*pp);
        
        z(ld+1:ld+vdim)=zz;
        ld=ld+vdim;

        zx(i,j,1)=zz(1,1);
        zy(i,j,1)=zz(2,1);
        zx(i,j,2)=zz(3,1);
        zy(i,j,2)=zz(4,1);
    end 

    % East nodes
    i=N+1;

    vdim=3;
    pdim=2;   
    pp=zeros(pdim,1);
    
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
    E_denom_E1_zz=4*JE1r3;
    JE2r2=abs(x5*(-y3 + y4) + x4*(y3 - y5) + x3*(-y4 + y5));
    E_denom_E2_zz=4*JE2r2;
    
    f=zeros(vdim,2);
    f(1,1)=-(1/2);
%     f(2,1)=0;
    f(3,1)=-(1/2);
%     f(1,2)=0;
    f(2,2)=-(1/2);
    f(3,2)=1/2;

    t=zeros(vdim,vdim);
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
    
    % CC.D en el nodo Este, fórmula de cuadratura para el término Pg de p en:
        % frontera Este 1
            Pgp_L7=int_dir_simpson(i,j,i,j-1,tt,3);
        % frontera Este 2
            Pgp_L14=int_dir_simpson(i,j,i,j+1,tt,3);
        Pgp=-(1/2)*[Pgp_L7;Pgp_L14;0];

    pp(1,1)=p(i-1,j-1);
    pp(2,1)=p(i-1,j);
    zz=t\(Pgp-f*pp);
    
    z(ld+1:ld+vdim)=zz;
    ld=ld+vdim;

    zx(i,j,1)=zz(1,1);
    zx(i,j,2)=zz(2,1);
    zy(i,j,2)=zz(3,1);
    
end

% North-West corner node
i=1;
j=N+1;
vdim=2;
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

Jer4=abs(x3*y1 - x4*y1 - x1*y3 + x4*y3 + x1*y4 - x3*y4);
NW_denom_zz=4*Jer4;

f=zeros(vdim,1);
f(1,1)=1/2;
f(2,1)=-(1/2);

t=zeros(vdim,vdim);
t(1,1)=Kinv(1,1)*(x3 - x4)^2 + (2*Kinv(1,2)*(x3 - x4) + Kinv(2,2)*(y3 - y4))*(y3 - y4);
t(2,1)=Kinv(1,1)*(x1 - x4)*(-x3 + x4) + Kinv(2,2)*(y1 - y4)*(-y3 + y4) + ...
       Kinv(1,2)*(x4*(y1 + y3 - 2*y4) + x3*(-y1 + y4) + x1*(-y3 + y4));
t(1,2)=t(2,1);
t(2,2)=Kinv(1,1)*(x1 - x4)^2 + (2*Kinv(1,2)*(x1 - x4) + Kinv(2,2)*(y1 - y4))*(y1 - y4);
t=t/NW_denom_zz;

% CC.D en el nodo N-W, fórmula de cuadratura para el término Pg de p en:
    % frontera Oeste
        Pgp_L18=int_dir_simpson(i,j,i,j-1,tt,3);
    % frontera Norte
        Pgp_L22=int_dir_simpson(i,j,i+1,j,tt,3);
    Pgp=(1/2)*[Pgp_L18;-Pgp_L22];

pp(1,1)=p(i,j-1);
zz=t\(Pgp-f*pp);

z(ld+1:ld+vdim)=zz;
ld=ld+vdim;

zx(i,j,1)=zz(1,1);
zy(i,j,1)=zz(2,1);

% North nodes (j=ny+1)
vdim=3;
pdim=2;
pp=zeros(pdim,1);

% Matrix f is the same for every North node
f=zeros(vdim,2);
f(1,1)=-(1/2);
% f(2,1)=0;
f(3,1)=-(1/2);
f(1,2)=1/2;
f(2,2)=-(1/2);
% f(3,2)=0;

t=zeros(vdim,vdim);

for i=2:N
    
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
    N_denom_E1_zz=4*JE1r3;
    JE2r4=abs(x6*(-y2 + y3) + x3*(y2 - y6) + x2*(-y3 + y6));
    N_denom_E2_zz=4*JE2r4;

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
    
    % CC.D en el nodo Norte, fórmula de cuadratura para el término Pg de p en:
        % frontera Norte 1
            Pgp_L23=int_dir_simpson(i,j,i+1,j,tt,3);
        % frontera Norte 2
            Pgp_L22=int_dir_simpson(i-1,j,i,j,tt,3);
        Pgp=-(1/2)*[0;Pgp_L23;Pgp_L22];
  
    pp(1,1)=p(i-1,j-1);
    pp(2,1)=p(i,j-1);
    zz=t\(Pgp-f*pp);
    
    z(ld+1:ld+vdim)=zz;
    ld=ld+vdim;

    zx(i,j,1)=zz(1,1);
    zy(i,j,1)=zz(2,1);
    zy(i,j,2)=zz(3,1);
end 

% North-East corner node (j=ny+1)
i=N+1;
vdim=2;
pdim=1;
pp=zeros(pdim,1);

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
NE_denom_zz=4*Jer3;

f=zeros(vdim,1);
f(1,1)=-(1/2);
f(2,1)=-(1/2);

t=zeros(vdim,vdim);
t(1,1)=Kinv(1,1)*(x3 - x4)^2 + (2*Kinv(1,2)*(x3 - x4) + Kinv(2,2)*(y3 - y4))*(y3 - y4);
t(2,1)=-(Kinv(1,1)*(x2 - x3)*(x3 - x4)) + Kinv(2,2)*(y2 - y3)*(y3 - y4) + ...
         Kinv(1,2)*(x4*(-y2 + y3) + x2*(y3 - y4) + x3*(y2 - 2*y3 + y4));
t(1,2)=t(2,1);
t(2,2)=Kinv(1,1)*(x2 - x3)^2 + (2*Kinv(1,2)*(x2 - x3) + Kinv(2,2)*(y2 - y3))*(y2 - y3);
t=t/NE_denom_zz;

% CC.D en el nodo N-E, fórmula de cuadratura para el término Pg de p en:
    % frontera Este
        Pgp_L21=int_dir_simpson(i,j,i,j-1,tt,3);
    % frontera Norte
        Pgp_L24=int_dir_simpson(i-1,j,i,j,tt,3);
    Pgp=-(1/2)*[Pgp_L21;Pgp_L24];

pp(1,1)=p(i-1,j-1);
zz=t\(Pgp-f*pp);

z(ld+1:ld+vdim)=zz;

zx(i,j,1)=zz(1,1);
zy(i,j,2)=zz(2,1);

return
end
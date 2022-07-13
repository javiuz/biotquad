function gamma=compute_gamma(u,p)

global NN x y alpha lambda mu 

N=NN;

gamma=zeros((N+1)*(N+1),1);

coef_d=alpha/(8*(lambda+mu));
aux_denom=16*mu*(lambda + mu);

% South-West corner node
i=1;
j=1;
ind2u=(i+(j-1)*N)*2;
ind1u=ind2u-1;
ind1p=ind2u/2;
vdim=4;
% vdim2=vdim/2;
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
denom=Jer1*aux_denom;

% Matriz local (A sigma, sigma)
a=zeros(vdim,vdim);
% a(1,1)=(lambda + 2*mu)*(x1 - x4)^2 + 2*(lambda + mu)*(y1 - y4)^2;
% a(1,2)=lambda*(-x1 + x4)*(y1 - y4);
% a(1,3)=(lambda + 2*mu)*(x1 - x2)*(x1 - x4) + 2*(lambda + mu)*(y1 - y2)*(y1 - y4);
% a(1,4)=lambda*(-x1 + x4)*(y1 - y2);
% a(2,1)=a(1,2);
% a(2,2)=2*(lambda + mu)*(x1 - x4)^2 + (lambda + 2*mu)*(y1 - y4)^2;
% a(2,3)=lambda*(-x1 + x2)*(y1 - y4);
% a(2,4)=2*(lambda + mu)*(x1 - x2)*(x1 - x4) + (lambda + 2*mu)*(y1 - y2)*(y1 - y4);
a(1,1)=denom;
a(2,2)=denom;
a(3,1)=(lambda + 2*mu)*(x1 - x2)*(x1 - x4) + 2*(lambda + mu)*(y1 - y2)*(y1 - y4);
a(3,2)=lambda*(-x1 + x2)*(y1 - y4);
a(3,3)=(lambda + 2*mu)*(x1 - x2)^2 + 2*(lambda + mu)*(y1 - y2)^2;
a(3,4)=lambda*(-x1 + x2)*(y1 - y2);
a(4,1)=lambda*(-x1 + x4)*(y1 - y2);
a(4,2)=2*(lambda + mu)*(x1 - x2)*(x1 - x4) + (lambda + 2*mu)*(y1 - y2)*(y1 - y4);
a(4,3)=a(3,4);
a(4,4)=2*(lambda + mu)*(x1 - x2)^2 + (lambda + 2*mu)*(y1 - y2)^2;
a=(1/denom)*a;

% Matriz local (A sigma, u)^t
% b=1/2*[-1 0;0 -1;-1 0;0 -1];
b=1/2*[0 0;0 0;-1 0;0 -1];

% Matriz local (A sigma, gamma)^t
% c=1/4*[y4-y1 x1-x4 y2-y1 x1-x2]';
c=1/4*[0 0 y2-y1 x1-x2]';

% Matriz local (A sigma, p)^t
d=zeros(vdim,1);
% d(1,1)=(-x1 + x4);
% d(2,1)=(-y1 + y4);
d(3,1)=(-x1 + x2);
d(4,1)=(-y1 + y2);
d=coef_d*d;

% Homogeneous Dirichlet (W) and Neuman (S) B.C.
gamma(pos)=(c'*(a\c))\(-c'*(a\b)*u(ind1u:ind2u)-c'*(a\d)*p(ind1p));  
pos=pos+1;

% South nodes (j=1)
vdim=6;

% Matrix b is the same for every South node
% b=(1/2)*[0 0 -1 0;0 0 0 -1;1 0 -1 0;0 1 0 -1;-1 0 0 0;0 -1 0 0];
b=(1/2)*[0 0 0 0;0 0 0 0;1 0 -1 0;0 1 0 -1;0 0 0 0;0 0 0 0];

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
    denom1=aux_denom*JE1r2;
    JE2r1=abs(x5*(-y2 + y3) + x3*(y2 - y5) + x2*(-y3 + y5));
    denom2=aux_denom*JE2r1;
    
%     a(1,1)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/denom2;
%     a(1,2)=-((lambda*(x2 - x3)*(y2 - y3))/denom2);
%     a(1,3)=((lambda + 2*mu)*(x2 - x3)*(x2 - x5) + 2*(lambda + mu)*(y2 - y3)*(y2 - y5))/denom2;
%     a(1,4)=-((lambda*(x2 - x3)*(y2 - y5))/denom2);
% %     a(1,5)=0;
% %     a(1,6)=0;
%     a(2,1)=a(1,2);
%     a(2,2)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/denom2;
%     a(2,3)=-((lambda*(x2 - x5)*(y2 - y3))/denom2);
%     a(2,4)=(2*(lambda + mu)*(x2 - x3)*(x2 - x5) + (lambda + 2*mu)*(y2 - y3)*(y2 - y5))/denom2;
% %     a(2,5)=0;
% %     a(2,6)=0;
    a(1,1)=1;
    a(2,2)=1;
    a(3,1)=((lambda + 2*mu)*(x2 - x3)*(x2 - x5) + 2*(lambda + mu)*(y2 - y3)*(y2 - y5))/denom2;
    a(3,2)=-((lambda*(x2 - x5)*(y2 - y3))/denom2);
    a(3,3)=((lambda + 2*mu)*(x1 - x2)^2 + 2*(lambda + mu)*(y1 - y2)^2)/denom1 + ...
           ((lambda + 2*mu)*(x2 - x5)^2 + 2*(lambda + mu)*(y2 - y5)^2)/denom2;
    a(3,4)=-((lambda*(x1 - x2)*(y1 - y2))/denom1) - (lambda*(x2 - x5)*(y2 - y5))/denom2;
    a(3,5)=((lambda + 2*mu)*(x1 - x2)*(x2 - x3) + 2*(lambda + mu)*(y1 - y2)*(y2 - y3))/denom1;
    a(3,6)=-((lambda*(x1 - x2)*(y2 - y3))/denom1);
    a(4,1)=-((lambda*(x2 - x3)*(y2 - y5))/denom2);
    a(4,2)=(2*(lambda + mu)*(x2 - x3)*(x2 - x5) + (lambda + 2*mu)*(y2 - y3)*(y2 - y5))/denom2;
    a(4,3)=a(3,4);
    a(4,4)=(2*(lambda + mu)*(x1 - x2)^2 + (lambda + 2*mu)*(y1 - y2)^2)/denom1 + ...
           (2*(lambda + mu)*(x2 - x5)^2 + (lambda + 2*mu)*(y2 - y5)^2)/denom2;
    a(4,5)=-((lambda*(x2 - x3)*(y1 - y2))/denom1);
    a(4,6)=(2*(lambda + mu)*(x1 - x2)*(x2 - x3) + (lambda + 2*mu)*(y1 - y2)*(y2 - y3))/denom1;
    a(5,5)=1;
    a(6,6)=1;
% %     a(5,1)=0;
% %     a(5,2)=0;
%     a(5,3)=a(3,5);
%     a(5,4)=a(4,5);
%     a(5,5)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/denom1;
%     a(5,6)=-((lambda*(x2 - x3)*(y2 - y3))/denom1);
% %     a(6,1)=0;
% %     a(6,2)=0;
%     a(6,3)=a(3,6);
%     a(6,4)=a(4,6);
%     a(6,5)=a(5,6);
%     a(6,6)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/denom1;
    
%     c=(1/4)*[y3-y2 x2-x3 y5-y1 x1-x5 y3-y2 x2-x3]';
    c=(1/4)*[0 0 y5-y1 x1-x5 0 0]';
    
%     d(1,1)=0;
%     d(2,1)=0;
    d(3,1)=(-x1 + x2);
    d(4,1)=(-y1 + y2);
%     d(5,1)=(-x2 + x3);
%     d(6,1)=(-y2 + y3);
%     d(1,2)=(-x2 + x3);
%     d(2,2)=(-y2 + y3);
    d(3,2)=(-x2 + x5);
    d(4,2)=(-y2 + y5);
%     d(5,2)=0;
%     d(6,2)=0;
    
    d=coef_d*d;
    
    % Homogeneous Neuman B.C.
    gamma(pos)=(c'*(a\c))\(-c'*(a\b)*u(ind1u:ind4u)-c'*(a\d)*p(ind1p:ind2p)); 
    pos=pos+1;
end

% South-East corner node (j=1)

% Homogeneous Neuman B.C.
% gamma(pos)=(c'*(a\c))\(-c'*(a\b)*u(ind1u:ind2u)-c'*(a\d)*p(ind1p)); 
gamma(pos)=0;
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
    denom1=aux_denom*JE1r4;
    JE2r1=abs(x6*(-y3 + y4) + x4*(y3 - y6) + x3*(-y4 + y6));
    denom2=aux_denom*JE2r1;
    
    a=zeros(vdim,vdim);    
    a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/denom1;
    a(1,2)=-((lambda*(x3 - x4)*(y3 - y4))/denom1);
    a(1,3)=-(((lambda + 2*mu)*(x1 - x4)*(x3 - x4) + 2*(lambda + mu)*(y1 - y4)*(y3 - y4))/denom1);
    a(1,4)=(lambda*(x3 - x4)*(y1 - y4))/denom1;
%     a(1,5)=0;
%     a(1,6)=0;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/denom1;
    a(2,3)=(lambda*(x1 - x4)*(y3 - y4))/denom1;
    a(2,4)=-((2*(lambda + mu)*(x1 - x4)*(x3 - x4) + (lambda + 2*mu)*(y1 - y4)*(y3 - y4))/denom1);
%     a(2,5)=0;
%     a(2,6)=0;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=((lambda + 2*mu)*(x1 - x4)^2 + 2*(lambda + mu)*(y1 - y4)^2)/denom1 + ...
           ((lambda + 2*mu)*(x4 - x6)^2 + 2*(lambda + mu)*(y4 - y6)^2)/denom2;
    a(3,4)=-((lambda*(x1 - x4)*(y1 - y4))/denom1) - (lambda*(x4 - x6)*(y4 - y6))/denom2;
    a(3,5)=-(((lambda + 2*mu)*(x3 - x4)*(x4 - x6) + 2*(lambda + mu)*(y3 - y4)*(y4 - y6))/denom2);
    a(3,6)=(lambda*(x4 - x6)*(y3 - y4))/denom2;
    a(4,1)=a(1,4);
    a(4,2)=a(2,4);
    a(4,3)=a(3,4);
    a(4,4)=(2*(lambda + mu)*(x1 - x4)^2 + (lambda + 2*mu)*(y1 - y4)^2)/denom1 + ...
           (2*(lambda + mu)*(x4 - x6)^2 + (lambda + 2*mu)*(y4 - y6)^2)/denom2;
    a(4,5)=(lambda*(x3 - x4)*(y4 - y6))/denom2;
    a(4,6)=(-2*(lambda + mu)*(x3 - x4)*(x4 - x6) - (lambda + 2*mu)*(y3 - y4)*(y4 - y6))/denom2;
%     a(5,1)=0;
%     a(5,2)=0;
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/denom2;
    a(5,6)=-((lambda*(x3 - x4)*(y3 - y4))/denom2);
%     a(6,1)=0;
%     a(6,2)=0;
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/denom2;
    
    b=(1/2)*[-1 0 0 0;0 -1 0 0;1 0 -1 0;0 1 0 -1;0 0 -1 0;0 0 0 -1];
    
    c=(1/4)*[y3-y4 x4-x3 y6-y1 x1-x6 y3-y4 x4-x3]';
    
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
    d(5,2)=(-x4 + x3);
    d(6,2)=(-y4 + y3);
    
    d=coef_d*d;
    
    % Homogeneous Dir.BC
    gamma(pos)=(c'*(a\c))\(-c'*(a\b)*u([ind1u:ind2u,ind3u:ind4u],1)-c'*(a\d)*p([ind1p,ind2p],1));
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
    denom1=aux_denom*JE1r3;
    JE2r4=abs(x6*(-y2 + y3) + x3*(y2 - y6) + x2*(-y3 + y6));
    denom2=aux_denom*JE2r4;
    JE3r1=abs(x8*(y3 - y6) + x3*(y6 - y8) + x6*(-y3 + y8));
    denom3=aux_denom*JE3r1;
    JE4r2=abs(x8*(-y3 + y4) + x4*(y3 - y8) + x3*(-y4 + y8));
    denom4=aux_denom*JE4r2;
    
    a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/denom1 + ...
           ((lambda + 2*mu)*(x3 - x6)^2 + 2*(lambda + mu)*(y3 - y6)^2)/denom2;
    a(1,2)=-((lambda*(x3 - x4)*(y3 - y4))/denom1) - (lambda*(x3 - x6)*(y3 - y6))/denom2;
    a(1,3)=((lambda + 2*mu)*(x2 - x3)*(x3 - x6) + 2*(lambda + mu)*(y2 - y3)*(y3 - y6))/denom2;
    a(1,4)=-((lambda*(x3 - x6)*(y2 - y3))/denom2);
%     a(1,5)=0;
%     a(1,6)=0;
    a(1,7)=-(((lambda + 2*mu)*(x2 - x3)*(x3 - x4) + 2*(lambda + mu)*(y2 - y3)*(y3 - y4))/denom1);
    a(1,8)=(lambda*(x3 - x4)*(y2 - y3))/denom1;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/denom1 + ...
           (2*(lambda + mu)*(x3 - x6)^2 + (lambda + 2*mu)*(y3 - y6)^2)/denom2;
    a(2,3)=-((lambda*(x2 - x3)*(y3 - y6))/denom2);
    a(2,4)=(2*(lambda + mu)*(x2 - x3)*(x3 - x6) + (lambda + 2*mu)*(y2 - y3)*(y3 - y6))/denom2;
%     a(2,5)=0;
%     a(2,6)=0;
    a(2,7)=(lambda*(x2 - x3)*(y3 - y4))/denom1;
    a(2,8)=(-2*(lambda + mu)*(x2 - x3)*(x3 - x4) - (lambda + 2*mu)*(y2 - y3)*(y3 - y4))/denom1;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/denom2 + ...
           ((lambda + 2*mu)*(x3 - x8)^2 + 2*(lambda + mu)*(y3 - y8)^2)/denom3;
    a(3,4)=-((lambda*(x2 - x3)*(y2 - y3))/denom2) - (lambda*(x3 - x8)*(y3 - y8))/denom3;
    a(3,5)=((lambda + 2*mu)*(x3 - x6)*(x3 - x8) + 2*(lambda + mu)*(y3 - y6)*(y3 - y8))/denom3;
    a(3,6)=-((lambda*(x3 - x8)*(y3 - y6))/denom3);
%     a(3,7)=0;
%     a(3,8)=0;
    a(4,1)=a(1,4);
    a(4,2)=a(2,4);
    a(4,3)=a(3,4);
    a(4,4)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/denom2 + ...
           (2*(lambda + mu)*(x3 - x8)^2 + (lambda + 2*mu)*(y3 - y8)^2)/denom3;
    a(4,5)=-((lambda*(x3 - x6)*(y3 - y8))/denom3);
    a(4,6)=(2*(lambda + mu)*(x3 - x6)*(x3 - x8) + (lambda + 2*mu)*(y3 - y6)*(y3 - y8))/denom3;
%     a(4,7)=0;
%     a(4,8)=0;
%     a(5,1)=0;
%     a(5,2)=0;
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/denom4 + ...
           ((lambda + 2*mu)*(x3 - x6)^2 + 2*(lambda + mu)*(y3 - y6)^2)/denom3;
    a(5,6)=-((lambda*(x3 - x4)*(y3 - y4))/denom4) - (lambda*(x3 - x6)*(y3 - y6))/denom3;
    a(5,7)=-(((lambda + 2*mu)*(x3 - x4)*(x3 - x8) + 2*(lambda + mu)*(y3 - y4)*(y3 - y8))/denom4);
    a(5,8)=(lambda*(x3 - x4)*(y3 - y8))/denom4;
%     a(6,1)=0;
%     a(6,2)=0;
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/denom4 + ...
           (2*(lambda + mu)*(x3 - x6)^2 + (lambda + 2*mu)*(y3 - y6)^2)/denom3;
    a(6,7)=(lambda*(x3 - x8)*(y3 - y4))/denom4;
    a(6,8)=-((2*(lambda + mu)*(x3 - x4)*(x3 - x8) + (lambda + 2*mu)*(y3 - y4)*(y3 - y8))/denom4);
    a(7,1)=a(1,7);
    a(7,2)=a(2,7);
%     a(7,3)=0;
%     a(7,4)=0;
    a(7,5)=a(5,7);
    a(7,6)=a(6,7);
    a(7,7)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/denom1 + ...
           ((lambda + 2*mu)*(x3 - x8)^2 + 2*(lambda + mu)*(y3 - y8)^2)/denom4;
    a(7,8)=-((lambda*(x2 - x3)*(y2 - y3))/denom1) - (lambda*(x3 - x8)*(y3 - y8))/denom4;
    a(8,1)=a(1,8);
    a(8,2)=a(2,8);
%     a(8,3)=0;
%     a(8,4)=0;
    a(8,5)=a(5,8);
    a(8,6)=a(6,8);
    a(8,7)=a(7,8);
    a(8,8)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/denom1 + ...
           (2*(lambda + mu)*(x3 - x8)^2 + (lambda + 2*mu)*(y3 - y8)^2)/denom4;
    
    c=(1/4)*[y6-y4 x4-x6 y8-y2 x2-x8 y6-y4 x4-x6 y8-y2 x2-x8]';
    
    % ¡¡OJO con las posiciones 3 y 4 (van intercaladas <- va antes la 4 que 
    % la 3!!
    d(1,1)=(x3 - x4);
    d(2,1)=(y3 - y4);
%     d(3,1)=0;
%     d(4,1)=0;
%     d(5,1)=0;
%     d(6,1)=0;
    d(7,1)=(-x2 + x3);
    d(8,1)=(-y2 + y3);
    d(1,2)=(x6 - x3);
    d(2,2)=(y6 - y3);
    d(3,2)=(-x2 + x3);
    d(4,2)=(-y2 + y3);
%     d(5,2)=0;
%     d(6,2)=0;
%     d(7,2)=0;
%     d(8,2)=0;

%     ¡OJO!: ¡Posicones 3 y 4 van intercaladas! (P4 va antes que P3)

%     d(1,3)=0;
%     d(2,3)=0;
%     d(3,3)=0;
%     d(4,3)=0;
    d(5,3)=(-x4 + x3);
    d(6,3)=(-y4 + y3);
    d(7,3)=(-x3 + x8);
    d(8,3)=(-y3 + y8);
    
%     d(1,4)=0;
%     d(2,4)=0;
    d(3,4)=(-x3 + x8);
    d(4,4)=(-y3 + y8);
    d(5,4)=(-x3 + x6);
    d(6,4)=(-y3 + y6);
%     d(7,4)=0;
%     d(8,4)=0;
    
    d=coef_d*d;
    
    % NO Dir.BC
    gamma(pos)=(c'*(a\c))\(-c'*(a\b)*u([ind1u:ind4u,ind5u:ind8u],1)-c'*(a\d)*p([ind1p:ind2p,ind3p:ind4p],1)); 
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
    denom1=aux_denom*JE1r3;
    JE2r2=abs(x5*(-y3 + y4) + x4*(y3 - y5) + x3*(-y4 + y5));
    denom2=aux_denom*JE2r2;
    
%     a=zeros(vdim,vdim);
    a=eye(vdim);
    
%     a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/denom1;
%     a(1,2)=-((lambda*(x3 - x4)*(y3 - y4))/denom1);
% %     a(1,3)=0;
% %     a(1,4)=0;
%     a(1,5)=-(((lambda + 2*mu)*(x2 - x3)*(x3 - x4) + 2*(lambda + mu)*(y2 - y3)*(y3 - y4))/denom1);
%     a(1,6)=(lambda*(x3 - x4)*(y2 - y3))/denom1;
%     a(2,1)=a(1,2);
%     a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/denom1;
% %     a(2,3)=0;
% %     a(2,4)=0:
%     a(2,5)=(lambda*(x2 - x3)*(y3 - y4))/denom1;
%     a(2,6)=-((2*(lambda + mu)*(x2 - x3)*(x3 - x4) + (lambda + 2*mu)*(y2 - y3)*(y3 - y4))/denom1);
% %     a(3,1)=0;
% %     a(3,2)=0;
%     a(3,3)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/denom2;
%     a(3,4)=-((lambda*(x3 - x4)*(y3 - y4))/denom2);
%     a(3,5)=-(((lambda + 2*mu)*(x3 - x4)*(x3 - x5) + 2*(lambda + mu)*(y3 - y4)*(y3 - y5))/denom2);
%     a(3,6)=(lambda*(x3 - x4)*(y3 - y5))/denom2;
% %     a(4,1)=0;
% %     a(4,2)=0;
%     a(4,3)=a(3,4);
%     a(4,4)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/denom2;
%     a(4,5)=(lambda*(x3 - x5)*(y3 - y4))/denom2;
%     a(4,6)=-((2*(lambda + mu)*(x3 - x4)*(x3 - x5) + (lambda + 2*mu)*(y3 - y4)*(y3 - y5))/denom2);
    a(5,1)=-(((lambda + 2*mu)*(x2 - x3)*(x3 - x4) + 2*(lambda + mu)*(y2 - y3)*(y3 - y4))/denom1);
    a(5,2)=(lambda*(x2 - x3)*(y3 - y4))/denom1;
    a(5,3)=-(((lambda + 2*mu)*(x3 - x4)*(x3 - x5) + 2*(lambda + mu)*(y3 - y4)*(y3 - y5))/denom2);
    a(5,4)=(lambda*(x3 - x5)*(y3 - y4))/denom2;
    a(5,5)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/denom1 + ...
           ((lambda + 2*mu)*(x3 - x5)^2 + 2*(lambda + mu)*(y3 - y5)^2)/denom2;
    a(5,6)=-((lambda*(x2 - x3)*(y2 - y3))/denom1) - (lambda*(x3 - x5)*(y3 - y5))/denom2;
    a(6,1)=(lambda*(x3 - x4)*(y2 - y3))/denom1;
    a(6,2)=-((2*(lambda + mu)*(x2 - x3)*(x3 - x4) + (lambda + 2*mu)*(y2 - y3)*(y3 - y4))/denom1);
    a(6,3)=(lambda*(x3 - x4)*(y3 - y5))/denom2;
    a(6,4)=-((2*(lambda + mu)*(x3 - x4)*(x3 - x5) + (lambda + 2*mu)*(y3 - y4)*(y3 - y5))/denom2);
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/denom1 +...
           (2*(lambda + mu)*(x3 - x5)^2 + (lambda + 2*mu)*(y3 - y5)^2)/denom2;
    
%     b=(1/2)*[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1;1 0 -1 0;0 1 0 -1];
    b=(1/2)*[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;1 0 -1 0;0 1 0 -1];
    
%     c=(1/4)*[y3-y4 x4-x3 y3-y4 x4-x3 y5-y2 x2-x5]';
    c=(1/4)*[0 0 0 0 y5-y2 x2-x5]';
    
    d=zeros(vdim,2);
%     d(1,1)=(x3 - x4);
%     d(2,1)=(y3 - y4);
% %     d(3,1)=0;
% %     d(4,1)=0;
    d(5,1)=(-x2 + x3);
    d(6,1)=(-y2 + y3);
% %     d(1,2)=0;
% %     d(2,2)=0;
%     d(3,2)=(-x4 + x3);
%     d(4,2)=(-y4 + y3);
    d(5,2)=(-x3 + x5);
    d(6,2)=(-y3 + y5);
    
    d=coef_d*d;
    
    % Homogeneous Neuman B.C.
    gamma(pos)=(c'*(a\c))\(-c'*(a\b)*u([ind1u:ind2u,ind3u:ind4u],1)-c'*(a\d)*p([ind1p,ind2p],1)); 
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
denom=aux_denom*Jer4;

a=zeros(vdim,vdim);
a(1,1)=(lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2;
a(1,2)=lambda*(-x3 + x4)*(y3 - y4);
a(1,3)=-((lambda + 2*mu)*(x1 - x4)*(x3 - x4)) - 2*(lambda + mu)*(y1 - y4)*(y3 - y4);
a(1,4)=lambda*(x3 - x4)*(y1 - y4);
a(2,1)=a(1,2);
a(2,2)=2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2;
a(2,3)=lambda*(x1 - x4)*(y3 - y4);
a(2,4)=-2*(lambda + mu)*(x1 - x4)*(x3 - x4) - (lambda + 2*mu)*(y1 - y4)*(y3 - y4);
a(3,3)=denom;
a(4,4)=denom;
% a(3,1)=a(1,3);
% a(3,2)=a(2,3);
% a(3,3)=(lambda + 2*mu)*(x1 - x4)^2 + 2*(lambda + mu)*(y1 - y4)^2;
% a(3,4)=lambda*(-x1 + x4)*(y1 - y4);
% a(4,1)=a(1,4);
% a(4,2)=a(2,4);
% a(4,3)=a(3,4);
% a(4,4)=2*(lambda + mu)*(x1 - x4)^2 + (lambda + 2*mu)*(y1 - y4)^2;
a=(1/denom)*a;

% b=(1/2)*[-1 0;0 -1;1 0;0 1];
b=(1/2)*[-1 0;0 -1;0 0;0 0];

% c=(1/4)*[y3-y4 x4-x3 y4-y1 x1-x4]';
c=(1/4)*[y3-y4 x4-x3 0 0]';

d=zeros(vdim,1);
d(1,1)=(x3 - x4);
d(2,1)=(y3 - y4);
% d(3,1)=(-x1 + x4);
% d(4,1)=(-y1 + y4);
d=coef_d*d;

% H.Dir.B.C. in the West part of the N-W node, and N.H.Neu.B.C. in the North boundary
        mod_L22=sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);
    Pgu_N=-mod_L22*[0;0;1;1];
    
g_term=-(1/4)*(y4-y1 + x1-x4)*mod_L22;    
% Non-homogeneous Dir.BC
gamma(pos)=(c'*(a\c))\(g_term + c'*(a\Pgu_N)-c'*(a\b)*u(ind1u:ind2u)-c'*(a\d)*p(ind1p)); 
pos=pos+1;

% North nodes (j=N+1)
vdim=6;

% Matrix b is the same for every North node
% b=(1/2)*[1 0 -1 0;0 1 0 -1;0 0 1 0;0 0 0 1;1 0 0 0;0 1 0 0];
b=(1/2)*[1 0 -1 0;0 1 0 -1;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];

% a=zeros(vdim,vdim);
a=eye(vdim);
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
    denom1=aux_denom*JE1r3;
    JE2r4=abs(x6*(-y2 + y3) + x3*(y2 - y6) + x2*(-y3 + y6));
    denom2=aux_denom*JE2r4;
    
    a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/denom1 + ...
           ((lambda + 2*mu)*(x3 - x6)^2 + 2*(lambda + mu)*(y3 - y6)^2)/denom2;
    a(1,2)=-((lambda*(x3 - x4)*(y3 - y4))/denom1) - (lambda*(x3 - x6)*(y3 - y6))/denom2;
    a(1,3)=((lambda + 2*mu)*(x2 - x3)*(x3 - x6) + 2*(lambda + mu)*(y2 - y3)*(y3 - y6))/denom2;
    a(1,4)=-((lambda*(x3 - x6)*(y2 - y3))/denom2);
    a(1,5)=-(((lambda + 2*mu)*(x2 - x3)*(x3 - x4) + 2*(lambda + mu)*(y2 - y3)*(y3 - y4))/denom1);
    a(1,6)=(lambda*(x3 - x4)*(y2 - y3))/denom1;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/denom1 + ...
           (2*(lambda + mu)*(x3 - x6)^2 + (lambda + 2*mu)*(y3 - y6)^2)/denom2;
    a(2,3)=-((lambda*(x2 - x3)*(y3 - y6))/denom2);
    a(2,4)=(2*(lambda + mu)*(x2 - x3)*(x3 - x6) + (lambda + 2*mu)*(y2 - y3)*(y3 - y6))/denom2;
    a(2,5)=(lambda*(x2 - x3)*(y3 - y4))/denom1;
    a(2,6)=(-2*(lambda + mu)*(x2 - x3)*(x3 - x4) - (lambda + 2*mu)*(y2 - y3)*(y3 - y4))/denom1;
%     a(3,1)=a(1,3);
%     a(3,2)=a(2,3);
%     a(3,3)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/denom2;
%     a(3,4)=-((lambda*(x2 - x3)*(y2 - y3))/denom2);
% %     a(3,5)=0;
% %     a(3,6)=0;
%     a(4,1)= a(1,4);
%     a(4,2)= a(2,4);
%     a(4,3)= a(3,4);
%     a(4,4)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/denom2;
% %     a(4,5)=0;
% %     a(4,6)=0;
%     a(5,1)=a(1,5);
%     a(5,2)=a(2,5);
% %     a(5,3)=0;
% %     a(5,4)=0;
%     a(5,5)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/denom1;
%     a(5,6)=-((lambda*(x2 - x3)*(y2 - y3))/denom1);
%     a(6,1)=a(1,6);
%     a(6,2)=a(2,6);
% %     a(6,3)=0;
% %     a(6,4)=0;
%     a(6,5)=a(5,6);
%     a(6,6)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/denom1;
    
%     c=(1/4)*[y6-y4 x4-x6 y3-y2 x2-x3 y3-y2 x2-x3]';
    c=(1/4)*[y6-y4 x4-x6 0 0 0 0]';
    
    d(1,1)=(x3 - x4);
    d(2,1)=(y3 - y4);
% %     d(3,1)=0;
% %     d(4,1)=0;
%     d(5,1)=(-x2 + x3);
%     d(6,1)=(-y2 + y3);
    d(1,2)=(x6 - x3);
    d(2,2)=(y6 - y3);
%     d(3,2)=(-x2 + x3);
%     d(4,2)=(-y2 + y3);
% %     d(5,2)=0;
% %     d(6,2)=0;

    d=coef_d*d;
    
    % Non-homogeneous Neumann B.C. in the North Node:
            mod_L23=sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);
            mod_L22=sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);
        Pgu_N=-[0;0;mod_L23;mod_L23;mod_L22;mod_L22];
        
        g_term=-(1/4)*(y3-y2 + x2-x3)*(mod_L23 + mod_L22);
      
    % Non-homogeneous Neuman B.C.
    gamma(pos)=(c'*(a\c))\(g_term + c'*(a\Pgu_N)-c'*(a\b)*u(ind1u:ind4u)-c'*(a\d)*p(ind1p:ind2p));
    pos=pos+1;
end

% North-East corner node (j=N+1)

% Non-homogeneous Neumann B.C. and homogeneous Neumann B.C. in the NE node:
% gamma(pos)=(c'*(a\c))\(g_term + c'*(a\Pgu_N)-c'*(a\b)*u(ind1u:ind2u)-c'*(a\d)*p(ind1p)); 
gamma(pos)=0; 
return
end
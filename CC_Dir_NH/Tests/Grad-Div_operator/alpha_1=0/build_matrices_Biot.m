function [M1,M2,M3,AspT,App]=build_matrices_Biot(delta_t)

global NN x y c0 alpha c_lambda c_mu

N=NN;

M1=sparse(2*N*N,2*N*N);
M2=sparse(2*N*N,N*N);
M3=sparse(N*N,N*N);
AspT=sparse(8*N*(N+1),N*N);
App=sparse(N*N,N*N);

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
x3=x(i+1,j+1);
y3=y(i+1,j+1);
x4=x(i,j+1);
y4=y(i,j+1);

xx=(x1+x2+x3+x4)/4;
yy=(y1+y2+y3+y4)/4;
E=sin(5*pi*xx)*sin(5*pi*yy)+5;

lambda=E*c_lambda;
mu=E*c_mu;

coef_d=alpha/(8*(lambda+mu));
coef_p=c0+alpha^2/(lambda+mu);

Jer1=abs(x4*(y1 - y2) + x1*(y2 - y4) + x2*(-y1 + y4));
denom=16*Jer1*mu*(lambda + mu);

% Matriz local (A sigma, sigma)
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
a=(1/denom)*a;

% Matriz local (A sigma, u)^t
b=1/2*[-1 0;0 -1;-1 0;0 -1];

% Matriz local (A sigma, gamma)^t
c=1/4*[y4-y1 x1-x4 y2-y1 x1-x2]';

% Matriz local (A sigma, p)^t
d=zeros(vdim,1);
d(1,1)=(-x1 + x4);
d(2,1)=(-y1 + y4);
d(3,1)=(-x1 + x2);
d(4,1)=(-y1 + y2);
d=coef_d*d;

% Matriz local (A z, p)^T
f=(1/2)*[1;1];

% Matriz local (A z, z)^{-1}
t=zeros(vdim2,vdim2);
t(1,1)=kinv(2,2,i,j,i,j);
t(1,2)=kinv(2,1,i,j,i,j);
t(2,1)=t(1,2);
t(2,2)=kinv(1,1,i,j,i,j);
t=(1/4)*t;

% Matriz local (A p, p)
k=coef_p;

localm1=b'*(a\b)-((b'*(a\c))*((c'*(a\c))\(c'*(a\b))));
localm2=b'*(a\d)-((b'*(a\c))*((c'*(a\c))\(c'*(a\d))));
% localm3=k-d'*(a\d)+(delta_t*f')*(t\f)+((d'*(a\c))*((c'*(a\c))\(c'*(a\d))));
localm3=k + (delta_t*f')*(t\f);

M1(ind1u:ind2u,ind1u:ind2u)=M1(ind1u:ind2u,ind1u:ind2u)+localm1;
M2(ind1u:ind2u,ind1p)=M2(ind1u:ind2u,ind1p)+localm2;
M3(ind1p,ind1p)=M3(ind1p,ind1p)+localm3;

indr=1:vdim;
AspT(indr,ind1p)=d;
ld=vdim;
App(ind1p,ind1p)=k;

% South nodes (j=1)
vdim=6;
vdim2=vdim/2;

% Matrices b and f are the same for every South node
b=(1/2)*[0 0 -1 0;0 0 0 -1;1 0 -1 0;0 1 0 -1;-1 0 0 0;0 -1 0 0];
f=(1/2)*[0 1;-1 1;1 0];

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
    x4=x(i-1,j+1);
    y4=y(i-1,j+1);
    x5=x(i+1,j);
    y5=y(i+1,j);
    x6=x(i+1,j+1);
    y6=y(i+1,j+1);
    
    xx1=(x1+x2+x3+x4)/4;
    yy1=(y1+y2+y3+y4)/4;
    E1=sin(5*pi*xx1)*sin(5*pi*yy1)+5;
    
    lambda1=E1*c_lambda;
    mu1=E1*c_mu;
    
    coef1_d=alpha/(8*(lambda1+mu1));
%     coef1_p=c0+alpha^2/(lambda1+mu1);
    
    xx2=(x5+x2+x3+x6)/4;
    yy2=(y5+y2+y3+y6)/4;
    E2=sin(5*pi*xx2)*sin(5*pi*yy2)+5;
    
    lambda2=E2*c_lambda;
    mu2=E2*c_mu;
    
    coef2_d=alpha/(8*(lambda2+mu2));
    coef2_p=c0+alpha^2/(lambda2+mu2);
    
    JE1r2=abs(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3));
    denom1=(16.*mu1*(lambda1 + mu1))*JE1r2;
    JE2r1=abs(x5*(-y2 + y3) + x3*(y2 - y5) + x2*(-y3 + y5));
    denom2=(16.*mu2*(lambda2 + mu2))*JE2r1;
    
    a(1,1)=((lambda2 + 2*mu2)*(x2 - x3)^2 + 2*(lambda2 + mu2)*(y2 - y3)^2)/denom2;
    a(1,2)=-((lambda2*(x2 - x3)*(y2 - y3))/denom2);
    a(1,3)=((lambda2 + 2*mu2)*(x2 - x3)*(x2 - x5) + 2*(lambda2 + mu2)*(y2 - y3)*(y2 - y5))/denom2;
    a(1,4)=-((lambda2*(x2 - x3)*(y2 - y5))/denom2);
%     a(1,5)=0;
%     a(1,6)=0;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda2 + mu2)*(x2 - x3)^2 + (lambda2 + 2*mu2)*(y2 - y3)^2)/denom2;
    a(2,3)=-((lambda2*(x2 - x5)*(y2 - y3))/denom2);
    a(2,4)=(2*(lambda2 + mu2)*(x2 - x3)*(x2 - x5) + (lambda2 + 2*mu2)*(y2 - y3)*(y2 - y5))/denom2;
%     a(2,5)=0;
%     a(2,6)=0;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=((lambda1 + 2*mu1)*(x1 - x2)^2 + 2*(lambda1 + mu1)*(y1 - y2)^2)/denom1 + ...
           ((lambda2 + 2*mu2)*(x2 - x5)^2 + 2*(lambda2 + mu2)*(y2 - y5)^2)/denom2;
    a(3,4)=-((lambda1*(x1 - x2)*(y1 - y2))/denom1) - (lambda2*(x2 - x5)*(y2 - y5))/denom2;
    a(3,5)=((lambda1 + 2*mu1)*(x1 - x2)*(x2 - x3) + 2*(lambda1 + mu1)*(y1 - y2)*(y2 - y3))/denom1;
    a(3,6)=-((lambda1*(x1 - x2)*(y2 - y3))/denom1);
    a(4,1)=a(1,4);
    a(4,2)=a(2,4);
    a(4,3)=a(3,4);
    a(4,4)=(2*(lambda1 + mu1)*(x1 - x2)^2 + (lambda1 + 2*mu1)*(y1 - y2)^2)/denom1 + ...
           (2*(lambda2 + mu2)*(x2 - x5)^2 + (lambda2 + 2*mu2)*(y2 - y5)^2)/denom2;
    a(4,5)=-((lambda1*(x2 - x3)*(y1 - y2))/denom1);
    a(4,6)=(2*(lambda1 + mu1)*(x1 - x2)*(x2 - x3) + (lambda1 + 2*mu1)*(y1 - y2)*(y2 - y3))/denom1;
%     a(5,1)=0;
%     a(5,2)=0;
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=((lambda1 + 2*mu1)*(x2 - x3)^2 + 2*(lambda1 + mu1)*(y2 - y3)^2)/denom1;
    a(5,6)=-((lambda1*(x2 - x3)*(y2 - y3))/denom1);
%     a(6,1)=0;
%     a(6,2)=0;
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda1 + mu1)*(x2 - x3)^2 + (lambda1 + 2*mu1)*(y2 - y3)^2)/denom1;
    
    c=(1/4)*[y3-y2 x2-x3 y5-y1 x1-x5 y3-y2 x2-x3]';
    
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
    
    d(:,1)=coef1_d*d(:,1);
    d(:,2)=coef2_d*d(:,2);

    t(1,1)=kinv(2,2,i,j,i,j);
    t(1,2)=kinv(2,1,i,j,i,j);
    %t(1,3)=0
    t(2,1)=t(1,2);
    t(2,2)=kinv(1,1,i-1,j,i,j)+kinv(1,1,i,j,i,j);
    t(2,3)=kinv(1,2,i-1,j,i,j);
    %t(3,1)=0
    t(3,2)=t(2,3);
    t(3,3)=kinv(2,2,i-1,j,i,j);
    t=(1/4)*t;
    
%     k1=coef1_p;
    k2=coef2_p;
    k=[0 0;0 k2];

    localm1=b'*(a\b)-((b'*(a\c))*((c'*(a\c))\(c'*(a\b))));
    localm2=b'*(a\d)-((b'*(a\c))*((c'*(a\c))\(c'*(a\d))));
%     localm3=k-d'*(a\d)+(delta_t*f')*(t\f)+((d'*(a\c))*((c'*(a\c))\(c'*(a\d))));
    localm3=k + (delta_t*f')*(t\f);
    
    M1(ind1u:ind4u,ind1u:ind4u)=M1(ind1u:ind4u,ind1u:ind4u)+localm1;   
    M2(ind1u:ind4u,ind1p:ind2p)=M2(ind1u:ind4u,ind1p:ind2p)+localm2;
    M3(ind1p:ind2p,ind1p:ind2p)=M3(ind1p:ind2p,ind1p:ind2p)+localm3;
    
    indr=ld+1:ld+vdim;
    AspT(indr,ind1p:ind2p)=d;
    ld=ld+vdim;
    App(ind2p,ind2p)=k2;
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
x4=x(i-1,j+1);
y4=y(i-1,j+1);

xx=(x1+x2+x3+x4)/4;
yy=(y1+y2+y3+y4)/4;
E=sin(5*pi*xx)*sin(5*pi*yy)+5;

lambda=E*c_lambda;
mu=E*c_mu;

coef_d=alpha/(8*(lambda+mu));
% coef_p=c0+alpha^2/(lambda+mu);

Jer2=abs(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3));
denom=16*Jer2*mu*(lambda + mu);

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
a=(1/denom)*a;

b=(1/2)*[1 0;0 1;-1 0;0 -1];

c=(1/4)*[y2-y1 x1-x2 y3-y2 x2-x3]';

d=zeros(vdim,1);
d(1,1)=(-x1 + x2);
d(2,1)=(-y1 + y2);
d(3,1)=(-x2 + x3);
d(4,1)=(-y2 + y3);
d=coef_d*d;

f=(1/2)*[-1;1];

t=zeros(vdim2,vdim2);
t(1,1)=kinv(1,1,i-1,j,i,j);
t(1,2)=kinv(1,2,i-1,j,i,j);
t(2,1)=t(1,2);
t(2,2)=kinv(2,2,i-1,j,i,j);
t=(1/4)*t;

% k=coef_p;

localm1=b'*(a\b)-((b'*(a\c))*((c'*(a\c))\(c'*(a\b))));
localm2=b'*(a\d)-((b'*(a\c))*((c'*(a\c))\(c'*(a\d))));
% localm3=-d'*(a\d)+(delta_t*f')*(t\f)+((d'*(a\c))*((c'*(a\c))\(c'*(a\d))));
localm3=(delta_t*f')*(t\f);

M1(ind1u:ind2u,ind1u:ind2u)=M1(ind1u:ind2u,ind1u:ind2u)+localm1;
M2(ind1u:ind2u,ind1p)=M2(ind1u:ind2u,ind1p)+localm2;
M3(ind1p,ind1p)=M3(ind1p,ind1p)+localm3;

indr=ld+1:ld+vdim;
AspT(indr,ind1p)=d;
ld=ld+vdim;

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
    x2=x(i+1,j-1);
    y2=y(i+1,j-1);
    x3=x(i+1,j);
    y3=y(i+1,j);
    x4=x(i,j);
    y4=y(i,j);
    x5=x(i+1,j+1);
    y5=y(i+1,j+1);
    x6=x(i,j+1);
    y6=y(i,j+1);
    
    xx1=(x1+x2+x3+x4)/4;
    yy1=(y1+y2+y3+y4)/4;
    E1=sin(5*pi*xx1)*sin(5*pi*yy1)+5;
    
    lambda1=E1*c_lambda;
    mu1=E1*c_mu;
    
    coef1_d=alpha/(8*(lambda1+mu1));
%     coef1_p=c0+alpha^2/(lambda1+mu1);
    
    xx2=(x3+x4+x5+x6)/4;
    yy2=(y3+y4+y5+y6)/4;
    E2=sin(5*pi*xx2)*sin(5*pi*yy2)+5;
    
    lambda2=E2*c_lambda;
    mu2=E2*c_mu;
    
    coef2_d=alpha/(8*(lambda2+mu2));
    coef2_p=c0+alpha^2/(lambda2+mu2);
    
    JE1r4=abs(x4*(y1 - y3) + x1*(y3 - y4) + x3*(-y1 + y4));
    denom1=(16.*mu1*(lambda1 + mu1))*JE1r4;
    JE2r1=abs(x6*(-y3 + y4) + x4*(y3 - y6) + x3*(-y4 + y6));
    denom2=(16.*mu2*(lambda2 + mu2))*JE2r1;
    
    a=zeros(vdim,vdim);    
    a(1,1)=((lambda1 + 2*mu1)*(x3 - x4)^2 + 2*(lambda1 + mu1)*(y3 - y4)^2)/denom1;
    a(1,2)=-((lambda1*(x3 - x4)*(y3 - y4))/denom1);
    a(1,3)=-(((lambda1 + 2*mu1)*(x1 - x4)*(x3 - x4) + 2*(lambda1 + mu1)*(y1 - y4)*(y3 - y4))/denom1);
    a(1,4)=(lambda1*(x3 - x4)*(y1 - y4))/denom1;
%     a(1,5)=0;
%     a(1,6)=0;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda1 + mu1)*(x3 - x4)^2 + (lambda1 + 2*mu1)*(y3 - y4)^2)/denom1;
    a(2,3)=(lambda1*(x1 - x4)*(y3 - y4))/denom1;
    a(2,4)=-((2*(lambda1 + mu1)*(x1 - x4)*(x3 - x4) + (lambda1 + 2*mu1)*(y1 - y4)*(y3 - y4))/denom1);
%     a(2,5)=0;
%     a(2,6)=0;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=((lambda1 + 2*mu1)*(x1 - x4)^2 + 2*(lambda1 + mu1)*(y1 - y4)^2)/denom1 + ...
           ((lambda2 + 2*mu2)*(x4 - x6)^2 + 2*(lambda2 + mu2)*(y4 - y6)^2)/denom2;
    a(3,4)=-((lambda1*(x1 - x4)*(y1 - y4))/denom1) - (lambda2*(x4 - x6)*(y4 - y6))/denom2;
    a(3,5)=-(((lambda2 + 2*mu2)*(x3 - x4)*(x4 - x6) + 2*(lambda2 + mu2)*(y3 - y4)*(y4 - y6))/denom2);
    a(3,6)=(lambda2*(x4 - x6)*(y3 - y4))/denom2;
    a(4,1)=a(1,4);
    a(4,2)=a(2,4);
    a(4,3)=a(3,4);
    a(4,4)=(2*(lambda1 + mu1)*(x1 - x4)^2 + (lambda1 + 2*mu1)*(y1 - y4)^2)/denom1 + ...
           (2*(lambda2 + mu2)*(x4 - x6)^2 + (lambda2 + 2*mu2)*(y4 - y6)^2)/denom2;
    a(4,5)=(lambda2*(x3 - x4)*(y4 - y6))/denom2;
    a(4,6)=(-2*(lambda2 + mu2)*(x3 - x4)*(x4 - x6) - (lambda2 + 2*mu2)*(y3 - y4)*(y4 - y6))/denom2;
%     a(5,1)=0;
%     a(5,2)=0;
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=((lambda2 + 2*mu2)*(x3 - x4)^2 + 2*(lambda2 + mu2)*(y3 - y4)^2)/denom2;
    a(5,6)=-((lambda2*(x3 - x4)*(y3 - y4))/denom2);
%     a(6,1)=0;
%     a(6,2)=0;
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda2 + mu2)*(x3 - x4)^2 + (lambda2 + 2*mu2)*(y3 - y4)^2)/denom2;
    
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
    
    d(:,1)=coef1_d*d(:,1);
    d(:,2)=coef2_d*d(:,2);
    
    f=(1/2)*[1 0;-1 1;0 1];
    
    t=zeros(vdim2,vdim2);  
    t(1,1)=kinv(1,1,i,j-1,i,j);
    t(1,2)=kinv(1,2,i,j-1,i,j);
    %t(1,3)=0;
    t(2,1)=t(1,2);
    t(2,2)=kinv(2,2,i,j-1,i,j)+kinv(2,2,i,j,i,j);
    t(2,3)=kinv(2,1,i,j,i,j);
    %t(3,1)=0;
    t(3,2)=t(2,3);
    t(3,3)=kinv(1,1,i,j,i,j);
    t=(1/4)*t;
    
%     k1=coef1_p;
    k2=coef2_p;
    k=[0 0;0 k2];
    
    localm1=b'*(a\b)-((b'*(a\c))*((c'*(a\c))\(c'*(a\b))));
    localm2=b'*(a\d)-((b'*(a\c))*((c'*(a\c))\(c'*(a\d))));
%     localm3=k-d'*(a\d)+(delta_t*f')*(t\f)+((d'*(a\c))*((c'*(a\c))\(c'*(a\d))));
    localm3=k + (delta_t*f')*(t\f);
    
    M1([ind1u:ind2u,ind3u:ind4u],[ind1u:ind2u,ind3u:ind4u])=M1([ind1u:ind2u,ind3u:ind4u],[ind1u:ind2u,ind3u:ind4u])+localm1;    
    M2([ind1u:ind2u,ind3u:ind4u],[ind1p,ind2p])=M2([ind1u:ind2u,ind3u:ind4u],[ind1p,ind2p])+localm2;    
    M3([ind1p,ind2p],[ind1p,ind2p])=M3([ind1p,ind2p],[ind1p,ind2p])+localm3;
    
    indr=ld+1:ld+vdim;
    AspT(indr,[ind1p,ind2p])=d;
    ld=ld+vdim;
    App(ind2p,ind2p)=k2;
    
    % Central nodes 
    vdim=8;
    vdim2=vdim/2;

    a=zeros(vdim,vdim);
    d=zeros(vdim,4);
    t=zeros(vdim2,vdim2);
    % Matrices b and f are the same for every central node
    b=(1/2)*[1 0 -1 0 0 0 0 0;0 1 0 -1 0 0 0 0;0 0 1 0 0 0 -1 0;...
             0 0 0 1 0 0 0 -1;0 0 0 0 1 0 -1 0;0 0 0 0 0 1 0 -1;...
             1 0 0 0 -1 0 0 0;0 1 0 0 0 -1 0 0];
    f=(1/2)*[-1 1 0 0;0 -1 0 1;0 0 -1 1;-1 0 1 0];
         
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
        
    x1=x(i-1,j-1);
    y1=y(i-1,j-1);
    x2=x(i,j-1);
    y2=y(i,j-1);
    x3=x(i,j);
    y3=y(i,j);
    x4=x(i-1,j);
    y4=y(i-1,j);
    x5=x(i+1,j-1);
    y5=y(i+1,j-1);
    x6=x(i+1,j);
    y6=y(i+1,j);
    x7=x(i+1,j+1);
    y7=y(i+1,j+1);
    x8=x(i,j+1);
    y8=y(i,j+1);
    x9=x(i-1,j+1);
    y9=y(i-1,j+1);
    
    xx1=(x1+x2+x3+x4)/4;
    yy1=(y1+y2+y3+y4)/4;
    E1=sin(5*pi*xx1)*sin(5*pi*yy1)+5;
    
    lambda1=E1*c_lambda;
    mu1=E1*c_mu;
    
    coef1_d=alpha/(8*(lambda1+mu1));
%     coef1_p=c0+alpha^2/(lambda1+mu1);
    
    xx2=(x5+x2+x3+x6)/4;
    yy2=(y5+y2+y3+y6)/4;
    E2=sin(5*pi*xx2)*sin(5*pi*yy2)+5;
    
    lambda2=E2*c_lambda;
    mu2=E2*c_mu;
    
    coef2_d=alpha/(8*(lambda2+mu2));
%     coef2_p=c0+alpha^2/(lambda2+mu2);
    
    % ¡OJO: Estos términos son de la sig.presión a P4!
    xx3=(x7+x8+x3+x6)/4;
    yy3=(y7+y8+y3+y6)/4;
    E3=sin(5*pi*xx3)*sin(5*pi*yy3)+5;
    
    lambda3=E3*c_lambda;
    mu3=E3*c_mu;
    
    coef3_d=alpha/(8*(lambda3+mu3));
    coef3_p=c0+alpha^2/(lambda3+mu3);
    
    % ¡OJO: Estos términos van una posición antes que P3!
    xx4=(x9+x8+x3+x4)/4;
    yy4=(y9+y8+y3+y4)/4;
    E4=sin(5*pi*xx4)*sin(5*pi*yy4)+5;
    
    lambda4=E4*c_lambda;
    mu4=E4*c_mu;
    
    coef4_d=alpha/(8*(lambda4+mu4));
%     coef4_p=c0+alpha^2/(lambda4+mu4);
    
    JE1r3=abs(x4*(y2 - y3) + x2*(y3 - y4) + x3*(-y2 + y4));
    denom1=(16.*mu1*(lambda1 + mu1))*JE1r3;
    JE2r4=abs(x6*(-y2 + y3) + x3*(y2 - y6) + x2*(-y3 + y6));
    denom2=(16.*mu2*(lambda2 + mu2))*JE2r4;
    JE3r1=abs(x8*(y3 - y6) + x3*(y6 - y8) + x6*(-y3 + y8));
    denom3=(16.*mu3*(lambda3 + mu3))*JE3r1;
    JE4r2=abs(x8*(-y3 + y4) + x4*(y3 - y8) + x3*(-y4 + y8));
    denom4=(16.*mu4*(lambda4 + mu4))*JE4r2;
    
    a(1,1)=((lambda1 + 2*mu1)*(x3 - x4)^2 + 2*(lambda1 + mu1)*(y3 - y4)^2)/denom1 + ...
           ((lambda2 + 2*mu2)*(x3 - x6)^2 + 2*(lambda2 + mu2)*(y3 - y6)^2)/denom2;
    a(1,2)=-((lambda1*(x3 - x4)*(y3 - y4))/denom1) - (lambda2*(x3 - x6)*(y3 - y6))/denom2;
    a(1,3)=((lambda2 + 2*mu2)*(x2 - x3)*(x3 - x6) + 2*(lambda2 + mu2)*(y2 - y3)*(y3 - y6))/denom2;
    a(1,4)=-((lambda2*(x3 - x6)*(y2 - y3))/denom2);
%     a(1,5)=0;
%     a(1,6)=0;
    a(1,7)=-(((lambda1 + 2*mu1)*(x2 - x3)*(x3 - x4) + 2*(lambda1 + mu1)*(y2 - y3)*(y3 - y4))/denom1);
    a(1,8)=(lambda1*(x3 - x4)*(y2 - y3))/denom1;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda1 + mu1)*(x3 - x4)^2 + (lambda1 + 2*mu1)*(y3 - y4)^2)/denom1 + ...
           (2*(lambda2 + mu2)*(x3 - x6)^2 + (lambda2 + 2*mu2)*(y3 - y6)^2)/denom2;
    a(2,3)=-((lambda2*(x2 - x3)*(y3 - y6))/denom2);
    a(2,4)=(2*(lambda2 + mu2)*(x2 - x3)*(x3 - x6) + (lambda2 + 2*mu2)*(y2 - y3)*(y3 - y6))/denom2;
%     a(2,5)=0;
%     a(2,6)=0;
    a(2,7)=(lambda1*(x2 - x3)*(y3 - y4))/denom1;
    a(2,8)=(-2*(lambda1 + mu1)*(x2 - x3)*(x3 - x4) - (lambda1 + 2*mu1)*(y2 - y3)*(y3 - y4))/denom1;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=((lambda2 + 2*mu2)*(x2 - x3)^2 + 2*(lambda2 + mu2)*(y2 - y3)^2)/denom2 + ...
           ((lambda3 + 2*mu3)*(x3 - x8)^2 + 2*(lambda3 + mu3)*(y3 - y8)^2)/denom3;
    a(3,4)=-((lambda2*(x2 - x3)*(y2 - y3))/denom2) - (lambda3*(x3 - x8)*(y3 - y8))/denom3;
    a(3,5)=((lambda3 + 2*mu3)*(x3 - x6)*(x3 - x8) + 2*(lambda3 + mu3)*(y3 - y6)*(y3 - y8))/denom3;
    a(3,6)=-((lambda3*(x3 - x8)*(y3 - y6))/denom3);
%     a(3,7)=0;
%     a(3,8)=0;
    a(4,1)=a(1,4);
    a(4,2)=a(2,4);
    a(4,3)=a(3,4);
    a(4,4)=(2*(lambda2 + mu2)*(x2 - x3)^2 + (lambda2 + 2*mu2)*(y2 - y3)^2)/denom2 + ...
           (2*(lambda3 + mu3)*(x3 - x8)^2 + (lambda3 + 2*mu3)*(y3 - y8)^2)/denom3;
    a(4,5)=-((lambda3*(x3 - x6)*(y3 - y8))/denom3);
    a(4,6)=(2*(lambda3 + mu3)*(x3 - x6)*(x3 - x8) + (lambda3 + 2*mu3)*(y3 - y6)*(y3 - y8))/denom3;
%     a(4,7)=0;
%     a(4,8)=0;
%     a(5,1)=0;
%     a(5,2)=0;
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=((lambda4 + 2*mu4)*(x3 - x4)^2 + 2*(lambda4 + mu4)*(y3 - y4)^2)/denom4 + ...
           ((lambda3 + 2*mu3)*(x3 - x6)^2 + 2*(lambda3 + mu3)*(y3 - y6)^2)/denom3;
    a(5,6)=-((lambda4*(x3 - x4)*(y3 - y4))/denom4) - (lambda3*(x3 - x6)*(y3 - y6))/denom3;
    a(5,7)=-(((lambda4 + 2*mu4)*(x3 - x4)*(x3 - x8) + 2*(lambda4 + mu4)*(y3 - y4)*(y3 - y8))/denom4);
    a(5,8)=(lambda4*(x3 - x4)*(y3 - y8))/denom4;
%     a(6,1)=0;
%     a(6,2)=0;
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda4 + mu4)*(x3 - x4)^2 + (lambda4 + 2*mu4)*(y3 - y4)^2)/denom4 + ...
           (2*(lambda3 + mu3)*(x3 - x6)^2 + (lambda3 + 2*mu3)*(y3 - y6)^2)/denom3;
    a(6,7)=(lambda4*(x3 - x8)*(y3 - y4))/denom4;
    a(6,8)=-((2*(lambda4 + mu4)*(x3 - x4)*(x3 - x8) + (lambda4 + 2*mu4)*(y3 - y4)*(y3 - y8))/denom4);
    a(7,1)=a(1,7);
    a(7,2)=a(2,7);
%     a(7,3)=0;
%     a(7,4)=0;
    a(7,5)=a(5,7);
    a(7,6)=a(6,7);
    a(7,7)=((lambda1 + 2*mu1)*(x2 - x3)^2 + 2*(lambda1 + mu1)*(y2 - y3)^2)/denom1 + ...
           ((lambda4 + 2*mu4)*(x3 - x8)^2 + 2*(lambda4 + mu4)*(y3 - y8)^2)/denom4;
    a(7,8)=-((lambda1*(x2 - x3)*(y2 - y3))/denom1) - (lambda4*(x3 - x8)*(y3 - y8))/denom4;
    a(8,1)=a(1,8);
    a(8,2)=a(2,8);
%     a(8,3)=0;
%     a(8,4)=0;
    a(8,5)=a(5,8);
    a(8,6)=a(6,8);
    a(8,7)=a(7,8);
    a(8,8)=(2*(lambda1 + mu1)*(x2 - x3)^2 + (lambda1 + 2*mu1)*(y2 - y3)^2)/denom1 + ...
           (2*(lambda4 + mu4)*(x3 - x8)^2 + (lambda4 + 2*mu4)*(y3 - y8)^2)/denom4;
    
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
    
    d(:,1)=coef1_d*d(:,1);
    d(:,2)=coef2_d*d(:,2);    
%     ¡OJO!: Índices 3 y 4 van intercalados    
    d(:,3)=coef4_d*d(:,3);
    d(:,4)=coef3_d*d(:,4);
    
    t(1,1)=kinv(1,1,i-1,j-1,i,j)+kinv(1,1,i,j-1,i,j);
    t(1,2)=kinv(1,2,i,j-1,i,j);
    %t(1,3)=0;
    t(1,4)=kinv(1,2,i-1,j-1,i,j);
    t(2,1)=t(1,2);
    t(2,2)=kinv(2,2,i,j-1,i,j)+kinv(2,2,i,j,i,j);
    t(2,3)=kinv(2,1,i,j,i,j);
    %t(2,4)=0;
    %t(3,1)=0;
    t(3,2)=t(2,3);
    t(3,3)=kinv(1,1,i,j,i,j)+kinv(1,1,i-1,j,i,j);
    t(3,4)=kinv(1,2,i-1,j,i,j);
    t(4,1)=t(1,4);
    %t(4,2)=0;
    t(4,3)=t(3,4);
    t(4,4)=kinv(2,2,i-1,j,i,j)+kinv(2,2,i-1,j-1,i,j);
    t=(1/4)*t;
    
%     k1=coef1_p;
%     k2=coef2_p;
%     SIGUE EL ORDEN DE LAS PRESIONES EN EL SISTEMA: p4 va antes que p3
%     k4=coef4_p;
    k3=coef3_p;
    k=[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 k3];
    
    localm1=b'*(a\b)-((b'*(a\c))*((c'*(a\c))\(c'*(a\b))));
    localm2=b'*(a\d)-((b'*(a\c))*((c'*(a\c))\(c'*(a\d))));
%     localm3=k-d'*(a\d)+(delta_t*f')*(t\f)+((d'*(a\c))*((c'*(a\c))\(c'*(a\d))));
    localm3=k + (delta_t*f')*(t\f);
    
    M1([ind1u:ind4u,ind5u:ind8u],[ind1u:ind4u,ind5u:ind8u])=M1([ind1u:ind4u,ind5u:ind8u],[ind1u:ind4u,ind5u:ind8u])+localm1;    
    M2([ind1u:ind4u,ind5u:ind8u],[ind1p:ind2p,ind3p:ind4p])=M2([ind1u:ind4u,ind5u:ind8u],[ind1p:ind2p,ind3p:ind4p])+localm2;    
    M3([ind1p:ind2p,ind3p:ind4p],[ind1p:ind2p,ind3p:ind4p])=M3([ind1p:ind2p,ind3p:ind4p],[ind1p:ind2p,ind3p:ind4p])+localm3;
    
    indr=ld+1:ld+vdim;
    AspT(indr,[ind1p:ind2p,ind3p:ind4p])=d;
    ld=ld+vdim;
    App(ind4p,ind4p)=k3;
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
    vdim2=vdim/2;
    
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
    x6=x(i-1,j+1);
    y6=y(i-1,j+1);   
    
    xx1=(x1+x2+x3+x4)/4;
    yy1=(y1+y2+y3+y4)/4;
    E1=sin(5*pi*xx1)*sin(5*pi*yy1)+5;
    
    lambda1=E1*c_lambda;
    mu1=E1*c_mu;
    
    coef1_d=alpha/(8*(lambda1+mu1));
%     coef1_p=c0+alpha^2/(lambda1+mu1);
    
    xx2=(x3+x4+x5+x6)/4;
    yy2=(y3+y4+y5+y6)/4;
    E2=sin(5*pi*xx2)*sin(5*pi*yy2)+5;
    
    lambda2=E2*c_lambda;
    mu2=E2*c_mu;
    
    coef2_d=alpha/(8*(lambda2+mu2));
%     coef2_p=c0+alpha^2/(lambda2+mu2);
    
    JE1r3=abs(x4*(y2 - y3) + x2*(y3 - y4) + x3*(-y2 + y4));
    denom1=(16.*mu1*(lambda1 + mu1))*JE1r3;
    JE2r2=abs(x5*(-y3 + y4) + x4*(y3 - y5) + x3*(-y4 + y5));
    denom2=(16.*mu2*(lambda2 + mu2))*JE2r2;
    
    a=zeros(vdim,vdim);
    a(1,1)=((lambda1 + 2*mu1)*(x3 - x4)^2 + 2*(lambda1 + mu1)*(y3 - y4)^2)/denom1;
    a(1,2)=-((lambda1*(x3 - x4)*(y3 - y4))/denom1);
%     a(1,3)=0;
%     a(1,4)=0;
    a(1,5)=-(((lambda1 + 2*mu1)*(x2 - x3)*(x3 - x4) + 2*(lambda1 + mu1)*(y2 - y3)*(y3 - y4))/denom1);
    a(1,6)=(lambda1*(x3 - x4)*(y2 - y3))/denom1;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda1 + mu1)*(x3 - x4)^2 + (lambda1 + 2*mu1)*(y3 - y4)^2)/denom1;
%     a(2,3)=0;
%     a(2,4)=0:
    a(2,5)=(lambda1*(x2 - x3)*(y3 - y4))/denom1;
    a(2,6)=-((2*(lambda1 + mu1)*(x2 - x3)*(x3 - x4) + (lambda1 + 2*mu1)*(y2 - y3)*(y3 - y4))/denom1);
%     a(3,1)=0;
%     a(3,2)=0;
    a(3,3)=((lambda2 + 2*mu2)*(x3 - x4)^2 + 2*(lambda2 + mu2)*(y3 - y4)^2)/denom2;
    a(3,4)=-((lambda2*(x3 - x4)*(y3 - y4))/denom2);
    a(3,5)=-(((lambda2 + 2*mu2)*(x3 - x4)*(x3 - x5) + 2*(lambda2 + mu2)*(y3 - y4)*(y3 - y5))/denom2);
    a(3,6)=(lambda2*(x3 - x4)*(y3 - y5))/denom2;
%     a(4,1)=0;
%     a(4,2)=0;
    a(4,3)=a(3,4);
    a(4,4)=(2*(lambda2 + mu2)*(x3 - x4)^2 + (lambda2 + 2*mu2)*(y3 - y4)^2)/denom2;
    a(4,5)=(lambda2*(x3 - x5)*(y3 - y4))/denom2;
    a(4,6)=-((2*(lambda2 + mu2)*(x3 - x4)*(x3 - x5) + (lambda2 + 2*mu2)*(y3 - y4)*(y3 - y5))/denom2);
    a(5,1)=a(1,5);
    a(5,2)=a(2,5);
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=((lambda1 + 2*mu1)*(x2 - x3)^2 + 2*(lambda1 + mu1)*(y2 - y3)^2)/denom1 + ...
           ((lambda2 + 2*mu2)*(x3 - x5)^2 + 2*(lambda2 + mu2)*(y3 - y5)^2)/denom2;
    a(5,6)=-((lambda1*(x2 - x3)*(y2 - y3))/denom1) - (lambda2*(x3 - x5)*(y3 - y5))/denom2;
    a(6,1)=a(1,6);
    a(6,2)=a(2,6);
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda1 + mu1)*(x2 - x3)^2 + (lambda1 + 2*mu1)*(y2 - y3)^2)/denom1 +...
           (2*(lambda2 + mu2)*(x3 - x5)^2 + (lambda2 + 2*mu2)*(y3 - y5)^2)/denom2;
    
    b=(1/2)*[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1;1 0 -1 0;0 1 0 -1];
    
    c=(1/4)*[y3-y4 x4-x3 y3-y4 x4-x3 y5-y2 x2-x5]';
    
    d=zeros(vdim,2);
    d(1,1)=(x3 - x4);
    d(2,1)=(y3 - y4);
%     d(3,1)=0;
%     d(4,1)=0;
    d(5,1)=(-x2 + x3);
    d(6,1)=(-y2 + y3);
%     d(1,2)=0;
%     d(2,2)=0;
    d(3,2)=(-x4 + x3);
    d(4,2)=(-y4 + y3);
    d(5,2)=(-x3 + x5);
    d(6,2)=(-y3 + y5);
    
    d(:,1)=coef1_d*d(:,1);
    d(:,2)=coef2_d*d(:,2);
    
    f=(1/2)*[-1 0;0 -1;-1 1];
    
    t=zeros(vdim2,vdim2);
    t(1,1)=kinv(1,1,i-1,j-1,i,j);
    %t(1,2)=0;
    t(1,3)=kinv(1,2,i-1,j-1,i,j);
    %t(2,1)=0;
    t(2,2)=kinv(1,1,i-1,j,i,j);
    t(2,3)=kinv(1,2,i-1,j,i,j);
    t(3,1)=t(1,3);
    t(3,2)=t(2,3);
    t(3,3)=kinv(2,2,i-1,j-1,i,j)+kinv(2,2,i-1,j,i,j);
    t=(1/4)*t;
    
%     k1=coef1_p;
%     k2=coef2_p;
%     k=[k1 0;0 k2];
    
    localm1=b'*(a\b)-((b'*(a\c))*((c'*(a\c))\(c'*(a\b))));
    localm2=b'*(a\d)-((b'*(a\c))*((c'*(a\c))\(c'*(a\d))));
%     localm3=-d'*(a\d)+(delta_t*f')*(t\f)+((d'*(a\c))*((c'*(a\c))\(c'*(a\d))));
    localm3=(delta_t*f')*(t\f);
    
    M1([ind1u:ind2u,ind3u:ind4u],[ind1u:ind2u,ind3u:ind4u])=M1([ind1u:ind2u,ind3u:ind4u],[ind1u:ind2u,ind3u:ind4u])+localm1; 
    M2([ind1u:ind2u,ind3u:ind4u],[ind1p,ind2p])=M2([ind1u:ind2u,ind3u:ind4u],[ind1p,ind2p])+localm2;
    M3([ind1p,ind2p],[ind1p,ind2p])=M3([ind1p,ind2p],[ind1p,ind2p])+localm3;
    
    indr=ld+1:ld+vdim;
    AspT(indr,[ind1p,ind2p])=d;
    ld=ld+vdim;
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
x2=x(i+1,j-1);
y2=y(i+1,j-1);
x3=x(i+1,j);
y3=y(i+1,j);
x4=x(i,j);
y4=y(i,j);

xx=(x1+x2+x3+x4)/4;
yy=(y1+y2+y3+y4)/4;
E=sin(5*pi*xx)*sin(5*pi*yy)+5;

lambda=E*c_lambda;
mu=E*c_mu;

coef_d=alpha/(8*(lambda+mu));
% coef_p=c0+alpha^2/(lambda+mu);

Jer4=abs(x4*(y1 - y3) + x1*(y3 - y4) + x3*(-y1 + y4));
denom=16*Jer4*mu*(lambda + mu);

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
a=(1/denom)*a;

b=(1/2)*[-1 0;0 -1;1 0;0 1];

c=(1/4)*[y3-y4 x4-x3 y4-y1 x1-x4]';

d=zeros(vdim,1);
d(1,1)=(x3 - x4);
d(2,1)=(y3 - y4);
d(3,1)=(-x1 + x4);
d(4,1)=(-y1 + y4);
d=coef_d*d;

f=(1/2)*[1;-1];

t=zeros(vdim2,vdim2);
t(1,1)=kinv(1,1,i,j-1,i,j);
t(1,2)=kinv(1,2,i,j-1,i,j);
t(2,1)=t(1,2);
t(2,2)=kinv(2,2,i,j-1,i,j);
t=(1/4)*t;

% k=coef_p;

localm1=b'*(a\b)-((b'*(a\c))*((c'*(a\c))\(c'*(a\b))));
localm2=b'*(a\d)-((b'*(a\c))*((c'*(a\c))\(c'*(a\d))));
% localm3=-d'*(a\d)+(delta_t*f')*(t\f)+((d'*(a\c))*((c'*(a\c))\(c'*(a\d))));
localm3=(delta_t*f')*(t\f);

M1(ind1u:ind2u,ind1u:ind2u)=M1(ind1u:ind2u,ind1u:ind2u)+localm1;
M2(ind1u:ind2u,ind1p)=M2(ind1u:ind2u,ind1p)+localm2;
M3(ind1p,ind1p)=M3(ind1p,ind1p)+localm3;

indr=ld+1:ld+vdim;
AspT(indr,ind1p)=d;
ld=ld+vdim;

% North nodes (j=N+1)
vdim=6;
vdim2=vdim/2;

% Matrices b and f are the same for every North node
b=(1/2)*[1 0 -1 0;0 1 0 -1;0 0 1 0;0 0 0 1;1 0 0 0;0 1 0 0];
f=(1/2)*[-1 1;0 -1;-1 0];

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
    x5=x(i+1,j-1);
    y5=y(i+1,j-1);
    x6=x(i+1,j);
    y6=y(i+1,j);   
    
    xx1=(x1+x2+x3+x4)/4;
    yy1=(y1+y2+y3+y4)/4;
    E1=sin(5*pi*xx1)*sin(5*pi*yy1)+5;
    
    lambda1=E1*c_lambda;
    mu1=E1*c_mu;
    
    coef1_d=alpha/(8*(lambda1+mu1));
%     coef1_p=c0+alpha^2/(lambda1+mu1);
    
    xx2=(x5+x2+x3+x6)/4;
    yy2=(y5+y2+y3+y6)/4;
    E2=sin(5*pi*xx2)*sin(5*pi*yy2)+5;
    
    lambda2=E2*c_lambda;
    mu2=E2*c_mu;
    
    coef2_d=alpha/(8*(lambda2+mu2));
%     coef2_p=c0+alpha^2/(lambda2+mu2);
    
    JE1r3=abs(x4*(y2 - y3) + x2*(y3 - y4) + x3*(-y2 + y4));
    denom1=(16.*mu1*(lambda1 + mu1))*JE1r3;
    JE2r4=abs(x6*(-y2 + y3) + x3*(y2 - y6) + x2*(-y3 + y6));
    denom2=(16.*mu2*(lambda2 + mu2))*JE2r4;
    
    a(1,1)=((lambda1 + 2*mu1)*(x3 - x4)^2 + 2*(lambda1 + mu1)*(y3 - y4)^2)/denom1 + ...
           ((lambda2 + 2*mu2)*(x3 - x6)^2 + 2*(lambda2 + mu2)*(y3 - y6)^2)/denom2;
    a(1,2)=-((lambda1*(x3 - x4)*(y3 - y4))/denom1) - (lambda2*(x3 - x6)*(y3 - y6))/denom2;
    a(1,3)=((lambda2 + 2*mu2)*(x2 - x3)*(x3 - x6) + 2*(lambda2 + mu2)*(y2 - y3)*(y3 - y6))/denom2;
    a(1,4)=-((lambda2*(x3 - x6)*(y2 - y3))/denom2);
    a(1,5)=-(((lambda1 + 2*mu1)*(x2 - x3)*(x3 - x4) + 2*(lambda1 + mu1)*(y2 - y3)*(y3 - y4))/denom1);
    a(1,6)=(lambda1*(x3 - x4)*(y2 - y3))/denom1;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda1 + mu1)*(x3 - x4)^2 + (lambda1 + 2*mu1)*(y3 - y4)^2)/denom1 + ...
           (2*(lambda2 + mu2)*(x3 - x6)^2 + (lambda2 + 2*mu2)*(y3 - y6)^2)/denom2;
    a(2,3)=-((lambda2*(x2 - x3)*(y3 - y6))/denom2);
    a(2,4)=(2*(lambda2 + mu2)*(x2 - x3)*(x3 - x6) + (lambda2 + 2*mu2)*(y2 - y3)*(y3 - y6))/denom2;
    a(2,5)=(lambda1*(x2 - x3)*(y3 - y4))/denom1;
    a(2,6)=(-2*(lambda1 + mu1)*(x2 - x3)*(x3 - x4) - (lambda1 + 2*mu1)*(y2 - y3)*(y3 - y4))/denom1;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=((lambda2 + 2*mu2)*(x2 - x3)^2 + 2*(lambda2 + mu2)*(y2 - y3)^2)/denom2;
    a(3,4)=-((lambda2*(x2 - x3)*(y2 - y3))/denom2);
%     a(3,5)=0;
%     a(3,6)=0;
    a(4,1)= a(1,4);
    a(4,2)= a(2,4);
    a(4,3)= a(3,4);
    a(4,4)=(2*(lambda2 + mu2)*(x2 - x3)^2 + (lambda2 + 2*mu2)*(y2 - y3)^2)/denom2;
%     a(4,5)=0;
%     a(4,6)=0;
    a(5,1)=a(1,5);
    a(5,2)=a(2,5);
%     a(5,3)=0;
%     a(5,4)=0;
    a(5,5)=((lambda1 + 2*mu1)*(x2 - x3)^2 + 2*(lambda1 + mu1)*(y2 - y3)^2)/denom1;
    a(5,6)=-((lambda1*(x2 - x3)*(y2 - y3))/denom1);
    a(6,1)=a(1,6);
    a(6,2)=a(2,6);
%     a(6,3)=0;
%     a(6,4)=0;
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda1 + mu1)*(x2 - x3)^2 + (lambda1 + 2*mu1)*(y2 - y3)^2)/denom1;
    
    c=(1/4)*[y6-y4 x4-x6 y3-y2 x2-x3 y3-y2 x2-x3]';
    
    d(1,1)=(x3 - x4);
    d(2,1)=(y3 - y4);
%     d(3,1)=0;
%     d(4,1)=0;
    d(5,1)=(-x2 + x3);
    d(6,1)=(-y2 + y3);
    d(1,2)=(x6 - x3);
    d(2,2)=(y6 - y3);
    d(3,2)=(-x2 + x3);
    d(4,2)=(-y2 + y3);
%     d(5,2)=0;
%     d(6,2)=0;

    d(:,1)=coef1_d*d(:,1);
    d(:,2)=coef2_d*d(:,2);

    t(1,1)=kinv(1,1,i-1,j-1,i,j)+kinv(1,1,i,j-1,i,j);
    t(1,2)=kinv(1,2,i,j-1,i,j);
    t(1,3)=kinv(1,2,i-1,j-1,i,j);
    t(2,1)=t(1,2);
    t(2,2)=kinv(2,2,i,j-1,i,j);
    %t(2,3)=0;
    t(3,1)=t(1,3);
    %t(3,2)=0;
    t(3,3)=kinv(2,2,i-1,j-1,i,j);
    t=(1/4)*t;

%     k1=coef1_p;
%     k2=coef2_p;
%     k=[k1 0;0 k2];
    
    localm1=b'*(a\b)-((b'*(a\c))*((c'*(a\c))\(c'*(a\b))));
    localm2=b'*(a\d)-((b'*(a\c))*((c'*(a\c))\(c'*(a\d))));
%     localm3=-d'*(a\d)+(delta_t*f')*(t\f)+((d'*(a\c))*((c'*(a\c))\(c'*(a\d))));
    localm3=(delta_t*f')*(t\f);
    
    M1(ind1u:ind4u,ind1u:ind4u)=M1(ind1u:ind4u,ind1u:ind4u)+localm1;
    M2(ind1u:ind4u,ind1p:ind2p)=M2(ind1u:ind4u,ind1p:ind2p)+localm2;  
    M3(ind1p:ind2p,ind1p:ind2p)=M3(ind1p:ind2p,ind1p:ind2p)+localm3;
    
    indr=ld+1:ld+vdim;
    AspT(indr,ind1p:ind2p)=d;
    ld=ld+vdim;
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

xx=(x1+x2+x3+x4)/4;
yy=(y1+y2+y3+y4)/4;
E=sin(5*pi*xx)*sin(5*pi*yy)+5;

lambda=E*c_lambda;
mu=E*c_mu;

coef_d=alpha/(8*(lambda+mu));
% coef_p=c0+alpha^2/(lambda+mu);

Jer3=abs(x4*(y2 - y3) + x2*(y3 - y4) + x3*(-y2 + y4));
denom=16*Jer3*mu*(lambda + mu);

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
a=(1/denom)*a;

b=(1/2)*[1 0;0 1;1 0;0 1];

c=(1/4)*[y3-y4 x4-x3 y3-y2 x2-x3]';

d=zeros(vdim,1);
d(1,1)=(x3 - x4);
d(2,1)=(y3 - y4);
d(3,1)=(-x2 + x3);
d(4,1)=(-y2 + y3);
d=coef_d*d;

f=(1/2)*[-1;-1];

t=zeros(vdim2,vdim2);
t(1,1)=kinv(1,1,i-1,j-1,i,j);
t(1,2)=kinv(1,2,i-1,j-1,i,j);
t(2,1)=t(1,2);
t(2,2)=kinv(2,2,i-1,j-1,i,j);
t=(1/4)*t;

% k=coef_p;

localm1=b'*(a\b)-((b'*(a\c))*((c'*(a\c))\(c'*(a\b))));
localm2=b'*(a\d)-((b'*(a\c))*((c'*(a\c))\(c'*(a\d))));
% localm3=-d'*(a\d)+(delta_t*f')*(t\f)+((d'*(a\c))*((c'*(a\c))\(c'*(a\d))));
localm3=(delta_t*f')*(t\f);

M1(ind1u:ind2u,ind1u:ind2u)=M1(ind1u:ind2u,ind1u:ind2u)+localm1;
M2(ind1u:ind2u,ind1p)=M2(ind1u:ind2u,ind1p)+localm2;
M3(ind1p,ind1p)=M3(ind1p,ind1p)+localm3;

indr=ld+1:ld+vdim;
AspT(indr,ind1p)=d;
return
end
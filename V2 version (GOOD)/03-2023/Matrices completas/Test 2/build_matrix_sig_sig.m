function matrix=build_matrix_sig_sig

global NN x y lambda mu

N=NN;

matrix=sparse(8*N*(N+1),8*N*(N+1));

% South-West corner node
i=1;
j=1;

vdim=4;
a=zeros(vdim,vdim);

x1=x(i,j);  
y1=y(i,j);            
x2=x(i+1,j);
y2=y(i+1,j);
% x3=x(i+1,j+1);
% y3=y(i+1,j+1);
x4=x(i,j+1);
y4=y(i,j+1);

Jer1=abs(x4*(y1 - y2) + x1*(y2 - y4) + x2*(-y1 + y4));
denom=16*Jer1*mu*(lambda + mu);

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

ind=1:vdim;
matrix(ind,ind)=a;
%     tf = issymmetric(a)
%     d = eig(a);
%     isposdef = all(d > 0)
ld=vdim;


% South nodes (j=1)
vdim=6;
a=zeros(vdim,vdim);

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
    
    JE1r2=abs(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3));
    JE2r1=abs(x5*(-y2 + y3) + x3*(y2 - y5) + x2*(-y3 + y5));
    
    a(1,1)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/JE2r1;
    a(1,2)=-((lambda*(x2 - x3)*(y2 - y3))/JE2r1);
    a(1,3)=((lambda + 2*mu)*(x2 - x3)*(x2 - x5) + 2*(lambda + mu)*(y2 - y3)*(y2 - y5))/JE2r1;
    a(1,4)=-((lambda*(x2 - x3)*(y2 - y5))/JE2r1);
%     a(1,5)=0;
%     a(1,6)=0;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/JE2r1;
    a(2,3)=-((lambda*(x2 - x5)*(y2 - y3))/JE2r1);
    a(2,4)=(2*(lambda + mu)*(x2 - x3)*(x2 - x5) + (lambda + 2*mu)*(y2 - y3)*(y2 - y5))/JE2r1;
%     a(2,5)=0;
%     a(2,6)=0;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=((lambda + 2*mu)*(x1 - x2)^2 + 2*(lambda + mu)*(y1 - y2)^2)/JE1r2 + ...
           ((lambda + 2*mu)*(x2 - x5)^2 + 2*(lambda + mu)*(y2 - y5)^2)/JE2r1;
    a(3,4)=-((lambda*(x1 - x2)*(y1 - y2))/JE1r2) - (lambda*(x2 - x5)*(y2 - y5))/JE2r1;
    a(3,5)=((lambda + 2*mu)*(x1 - x2)*(x2 - x3) + 2*(lambda + mu)*(y1 - y2)*(y2 - y3))/JE1r2;
    a(3,6)=-((lambda*(x1 - x2)*(y2 - y3))/JE1r2);
    a(4,1)=a(1,4);
    a(4,2)=a(2,4);
    a(4,3)=a(3,4);
    a(4,4)=(2*(lambda + mu)*(x1 - x2)^2 + (lambda + 2*mu)*(y1 - y2)^2)/JE1r2 + ...
           (2*(lambda + mu)*(x2 - x5)^2 + (lambda + 2*mu)*(y2 - y5)^2)/JE2r1;
    a(4,5)=-((lambda*(x2 - x3)*(y1 - y2))/JE1r2);
    a(4,6)=(2*(lambda + mu)*(x1 - x2)*(x2 - x3) + (lambda + 2*mu)*(y1 - y2)*(y2 - y3))/JE1r2;
%     a(5,1)=0;
%     a(5,2)=0;
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/JE1r2;
    a(5,6)=-((lambda*(x2 - x3)*(y2 - y3))/JE1r2);
%     a(6,1)=0;
%     a(6,2)=0;
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/JE1r2;
    a=(1/(16.*mu*(lambda + mu)))*a;
    
    ind=ld+1:ld+vdim;
    matrix(ind,ind)=a;
%         tf = issymmetric(a)
%         d = eig(a);
%         isposdef = all(d > 0)
    ld=ld+vdim;
end

% South-East corner node (j=1)
i=N+1;

vdim=4;
a=zeros(vdim,vdim);

x1=x(i-1,j);  
y1=y(i-1,j);            
x2=x(i,j);
y2=y(i,j);
x3=x(i,j+1);
y3=y(i,j+1);
% x4=x(i-1,j+1);
% y4=y(i-1,j+1);

Jer2=abs(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3));
denom=16*Jer2*mu*(lambda + mu);

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

ind=ld+1:ld+vdim;
matrix(ind,ind)=a;
%     tf = issymmetric(a)
%     d = eig(a);
%     isposdef = all(d > 0)
ld=ld+vdim;


for j=2:N

    % West nodes
    i=1;
    
    vdim=6;
    a=zeros(vdim,vdim);  
   
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
    JE2r1=abs(x6*(-y3 + y4) + x4*(y3 - y6) + x3*(-y4 + y6));
    
    a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/JE1r4;
    a(1,2)=-((lambda*(x3 - x4)*(y3 - y4))/JE1r4);
    a(1,3)=-(((lambda + 2*mu)*(x1 - x4)*(x3 - x4) + 2*(lambda + mu)*(y1 - y4)*(y3 - y4))/JE1r4);
    a(1,4)=(lambda*(x3 - x4)*(y1 - y4))/JE1r4;
%     a(1,5)=0;
%     a(1,6)=0;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/JE1r4;
    a(2,3)=(lambda*(x1 - x4)*(y3 - y4))/JE1r4;
    a(2,4)=-((2*(lambda + mu)*(x1 - x4)*(x3 - x4) + (lambda + 2*mu)*(y1 - y4)*(y3 - y4))/JE1r4);
%     a(2,5)=0;
%     a(2,6)=0;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=((lambda + 2*mu)*(x1 - x4)^2 + 2*(lambda + mu)*(y1 - y4)^2)/JE1r4 + ...
           ((lambda + 2*mu)*(x4 - x6)^2 + 2*(lambda + mu)*(y4 - y6)^2)/JE2r1;
    a(3,4)=-((lambda*(x1 - x4)*(y1 - y4))/JE1r4) - (lambda*(x4 - x6)*(y4 - y6))/JE2r1;
    a(3,5)=-(((lambda + 2*mu)*(x3 - x4)*(x4 - x6) + 2*(lambda + mu)*(y3 - y4)*(y4 - y6))/JE2r1);
    a(3,6)=(lambda*(x4 - x6)*(y3 - y4))/JE2r1;
    a(4,1)=a(1,4);
    a(4,2)=a(2,4);
    a(4,3)=a(3,4);
    a(4,4)=(2*(lambda + mu)*(x1 - x4)^2 + (lambda + 2*mu)*(y1 - y4)^2)/JE1r4 + ...
           (2*(lambda + mu)*(x4 - x6)^2 + (lambda + 2*mu)*(y4 - y6)^2)/JE2r1;
    a(4,5)=(lambda*(x3 - x4)*(y4 - y6))/JE2r1;
    a(4,6)=(-2*(lambda + mu)*(x3 - x4)*(x4 - x6) - (lambda + 2*mu)*(y3 - y4)*(y4 - y6))/JE2r1;
%     a(5,1)=0;
%     a(5,2)=0;
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/JE2r1;
    a(5,6)=-((lambda*(x3 - x4)*(y3 - y4))/JE2r1);
%     a(6,1)=0;
%     a(6,2)=0;
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/JE2r1;
    a=(1/(16.*mu*(lambda + mu)))*a;
    
    ind=ld+1:ld+vdim;
    matrix(ind,ind)=a;
%         tf = issymmetric(a)
%         d = eig(a);
%         isposdef = all(d > 0)
    ld=ld+vdim;
    
     % Central nodes 
    vdim=8;

    a=zeros(vdim,vdim);
    
    for i=2:N
        
%     x1=x(i - 1,j - 1);
%     y1=y(i - 1,j - 1);
    x2=x(i, j - 1);
    y2=y(i, j - 1);
    x3=x(i, j);
    y3=y(i, j);
    x4=x(i - 1, j);
    y4=y(i - 1, j);
%     x5=x(i + 1, j - 1);
%     y5=y(i + 1, j - 1);
    x6=x(i + 1, j);
    y6=y(i + 1, j);
%     x7=x(i + 1, j + 1);
%     y7=y(i + 1, j + 1);
    x8=x(i, j + 1);
    y8=y(i, j + 1);
%     x9=x(i - 1, j + 1);
%     y9=y(i - 1, j + 1);  
    
    JE1r3=abs(x4*(y2 - y3) + x2*(y3 - y4) + x3*(-y2 + y4));
    JE2r4=abs(x6*(-y2 + y3) + x3*(y2 - y6) + x2*(-y3 + y6));
    JE3r1=abs(x8*(y3 - y6) + x3*(y6 - y8) + x6*(-y3 + y8));
    JE4r2=abs(x8*(-y3 + y4) + x4*(y3 - y8) + x3*(-y4 + y8));
    
    a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/JE1r3 + ...
           ((lambda + 2*mu)*(x3 - x6)^2 + 2*(lambda + mu)*(y3 - y6)^2)/JE2r4;
    a(1,2)=-((lambda*(x3 - x4)*(y3 - y4))/JE1r3) - (lambda*(x3 - x6)*(y3 - y6))/JE2r4;
    a(1,3)=((lambda + 2*mu)*(x2 - x3)*(x3 - x6) + 2*(lambda + mu)*(y2 - y3)*(y3 - y6))/JE2r4;
    a(1,4)=-((lambda*(x3 - x6)*(y2 - y3))/JE2r4);
%     a(1,5)=0;
%     a(1,6)=0;
    a(1,7)=-(((lambda + 2*mu)*(x2 - x3)*(x3 - x4) + 2*(lambda + mu)*(y2 - y3)*(y3 - y4))/JE1r3);
    a(1,8)=(lambda*(x3 - x4)*(y2 - y3))/JE1r3;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/JE1r3 + ...
           (2*(lambda + mu)*(x3 - x6)^2 + (lambda + 2*mu)*(y3 - y6)^2)/JE2r4;
    a(2,3)=-((lambda*(x2 - x3)*(y3 - y6))/JE2r4);
    a(2,4)=(2*(lambda + mu)*(x2 - x3)*(x3 - x6) + (lambda + 2*mu)*(y2 - y3)*(y3 - y6))/JE2r4;
%     a(2,5)=0;
%     a(2,6)=0;
    a(2,7)=(lambda*(x2 - x3)*(y3 - y4))/JE1r3;
    a(2,8)=(-2*(lambda + mu)*(x2 - x3)*(x3 - x4) - (lambda + 2*mu)*(y2 - y3)*(y3 - y4))/JE1r3;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/JE2r4 + ...
           ((lambda + 2*mu)*(x3 - x8)^2 + 2*(lambda + mu)*(y3 - y8)^2)/JE3r1;
    a(3,4)=-((lambda*(x2 - x3)*(y2 - y3))/JE2r4) - (lambda*(x3 - x8)*(y3 - y8))/JE3r1;
    a(3,5)=((lambda + 2*mu)*(x3 - x6)*(x3 - x8) + 2*(lambda + mu)*(y3 - y6)*(y3 - y8))/JE3r1;
    a(3,6)=-((lambda*(x3 - x8)*(y3 - y6))/JE3r1);
%     a(3,7)=0;
%     a(3,8)=0;
    a(4,1)=a(1,4);
    a(4,2)=a(2,4);
    a(4,3)=a(3,4);
    a(4,4)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/JE2r4 + ...
           (2*(lambda + mu)*(x3 - x8)^2 + (lambda + 2*mu)*(y3 - y8)^2)/JE3r1;
    a(4,5)=-((lambda*(x3 - x6)*(y3 - y8))/JE3r1);
    a(4,6)=(2*(lambda + mu)*(x3 - x6)*(x3 - x8) + (lambda + 2*mu)*(y3 - y6)*(y3 - y8))/JE3r1;
%     a(4,7)=0;
%     a(4,8)=0;
%     a(5,1)=0;
%     a(5,2)=0;
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/JE4r2 + ...
           ((lambda + 2*mu)*(x3 - x6)^2 + 2*(lambda + mu)*(y3 - y6)^2)/JE3r1;
    a(5,6)=-((lambda*(x3 - x4)*(y3 - y4))/JE4r2) - (lambda*(x3 - x6)*(y3 - y6))/JE3r1;
    a(5,7)=-(((lambda + 2*mu)*(x3 - x4)*(x3 - x8) + 2*(lambda + mu)*(y3 - y4)*(y3 - y8))/JE4r2);
    a(5,8)=(lambda*(x3 - x4)*(y3 - y8))/JE4r2;
%     a(6,1)=0;
%     a(6,2)=0;
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/JE4r2 + ...
           (2*(lambda + mu)*(x3 - x6)^2 + (lambda + 2*mu)*(y3 - y6)^2)/JE3r1;
    a(6,7)=(lambda*(x3 - x8)*(y3 - y4))/JE4r2;
    a(6,8)=-((2*(lambda + mu)*(x3 - x4)*(x3 - x8) + (lambda + 2*mu)*(y3 - y4)*(y3 - y8))/JE4r2);
    a(7,1)=a(1,7);
    a(7,2)=a(2,7);
%     a(7,3)=0;
%     a(7,4)=0;
    a(7,5)=a(5,7);
    a(7,6)=a(6,7);
    a(7,7)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/JE1r3 + ...
           ((lambda + 2*mu)*(x3 - x8)^2 + 2*(lambda + mu)*(y3 - y8)^2)/JE4r2;
    a(7,8)=-((lambda*(x2 - x3)*(y2 - y3))/JE1r3) - (lambda*(x3 - x8)*(y3 - y8))/JE4r2;
    a(8,1)=a(1,8);
    a(8,2)=a(2,8);
%     a(8,3)=0;
%     a(8,4)=0;
    a(8,5)=a(5,8);
    a(8,6)=a(6,8);
    a(8,7)=a(7,8);
    a(8,8)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/JE1r3 + ...
           (2*(lambda + mu)*(x3 - x8)^2 + (lambda + 2*mu)*(y3 - y8)^2)/JE4r2;
    a=(1/(16.*mu*(lambda + mu)))*a;
    
    ind=ld+1:ld+vdim;
    matrix(ind,ind)=a;
%         tf = issymmetric(a)
%         d = eig(a);
%         isposdef = all(d > 0)
    ld=ld+vdim;
    end
    
    % East nodes
    i=N+1;
  
    vdim=6;
    a=zeros(vdim,vdim);
    
%     x1=x(i - 1, j - 1);
%     y1=y(i - 1, j - 1);
    x2=x(i, j - 1);
    y2=y(i, j - 1);
    x3=x(i, j);
    y3=y(i, j);
    x4=x(i - 1, j);
    y4=y(i - 1, j);
    x5=x(i, j + 1);
    y5=y(i, j + 1);
%     x6=x(i - 1, j + 1);
%     y6=y(i - 1, j + 1);
    
    JE1r3=abs(x4*(y2 - y3) + x2*(y3 - y4) + x3*(-y2 + y4));
    JE2r2=abs(x5*(-y3 + y4) + x4*(y3 - y5) + x3*(-y4 + y5));
    
    a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/JE1r3;
    a(1,2)=-((lambda*(x3 - x4)*(y3 - y4))/JE1r3);
%     a(1,3)=0;
%     a(1,4)=0;
    a(1,5)=-(((lambda + 2*mu)*(x2 - x3)*(x3 - x4) + 2*(lambda + mu)*(y2 - y3)*(y3 - y4))/JE1r3);
    a(1,6)=(lambda*(x3 - x4)*(y2 - y3))/JE1r3;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/JE1r3;
%     a(2,3)=0;
%     a(2,4)=0:
    a(2,5)=(lambda*(x2 - x3)*(y3 - y4))/JE1r3;
    a(2,6)=-((2*(lambda + mu)*(x2 - x3)*(x3 - x4) + (lambda + 2*mu)*(y2 - y3)*(y3 - y4))/JE1r3);
%     a(3,1)=0;
%     a(3,2)=0;
    a(3,3)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/JE2r2;
    a(3,4)=-((lambda*(x3 - x4)*(y3 - y4))/JE2r2);
    a(3,5)=-(((lambda + 2*mu)*(x3 - x4)*(x3 - x5) + 2*(lambda + mu)*(y3 - y4)*(y3 - y5))/JE2r2);
    a(3,6)=(lambda*(x3 - x4)*(y3 - y5))/JE2r2;
%     a(4,1)=0;
%     a(4,2)=0;
    a(4,3)=a(3,4);
    a(4,4)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/JE2r2;
    a(4,5)=(lambda*(x3 - x5)*(y3 - y4))/JE2r2;
    a(4,6)=-((2*(lambda + mu)*(x3 - x4)*(x3 - x5) + (lambda + 2*mu)*(y3 - y4)*(y3 - y5))/JE2r2);
    a(5,1)=a(1,5);
    a(5,2)=a(2,5);
    a(5,3)=a(3,5);
    a(5,4)=a(4,5);
    a(5,5)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/JE1r3 + ...
           ((lambda + 2*mu)*(x3 - x5)^2 + 2*(lambda + mu)*(y3 - y5)^2)/JE2r2;
    a(5,6)=-((lambda*(x2 - x3)*(y2 - y3))/JE1r3) - (lambda*(x3 - x5)*(y3 - y5))/JE2r2;
    a(6,1)=a(1,6);
    a(6,2)=a(2,6);
    a(6,3)=a(3,6);
    a(6,4)=a(4,6);
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/JE1r3 +...
           (2*(lambda + mu)*(x3 - x5)^2 + (lambda + 2*mu)*(y3 - y5)^2)/JE2r2;
    a=(1/(16.*mu*(lambda + mu)))*a;
    
    ind=ld+1:ld+vdim;
    matrix(ind,ind)=a;
%          tf = issymmetric(a)
%          d = eig(a);
%          isposdef = all(d > 0)
    ld=ld+vdim;
end

% North-West corner node
i=1;
j=N+1;

vdim=4;
a=zeros(vdim,vdim);

x1=x(i, j - 1);
y1=y(i, j - 1);           
% x2=x(i + 1, j - 1);
% y2=y(i + 1, j - 1);
x3=x(i + 1, j);
y3=y(i + 1, j);
x4=x(i, j);
y4=y(i, j);

Jer4=abs(x4*(y1 - y3) + x1*(y3 - y4) + x3*(-y1 + y4));
denom=16*Jer4*mu*(lambda + mu);

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

ind=ld+1:ld+vdim;
matrix(ind,ind)=a;
%     tf = issymmetric(a)
%     d = eig(a);
%     isposdef = all(d > 0)
ld=ld+vdim;

% North nodes (j=N+1)
vdim=6;
a=zeros(vdim,vdim);

for i=2:N
    
%     x1=x(i - 1, j - 1);
%     y1=y(i - 1, j - 1);
    x2=x(i, j - 1);
    y2=y(i, j - 1);
    x3=x(i, j);
    y3=y(i, j);
    x4=x(i - 1, j);
    y4=y(i - 1, j);
%     x5=x(i + 1, j - 1);
%     y5=y(i + 1, j - 1);
    x6=x(i + 1, j);
    y6=y(i + 1, j);
    
    JE1r3=abs(x4*(y2 - y3) + x2*(y3 - y4) + x3*(-y2 + y4));
    JE2r4=abs(x6*(-y2 + y3) + x3*(y2 - y6) + x2*(-y3 + y6));
    
    a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2)/JE1r3 + ...
           ((lambda + 2*mu)*(x3 - x6)^2 + 2*(lambda + mu)*(y3 - y6)^2)/JE2r4;
    a(1,2)=-((lambda*(x3 - x4)*(y3 - y4))/JE1r3) - (lambda*(x3 - x6)*(y3 - y6))/JE2r4;
    a(1,3)=((lambda + 2*mu)*(x2 - x3)*(x3 - x6) + 2*(lambda + mu)*(y2 - y3)*(y3 - y6))/JE2r4;
    a(1,4)=-((lambda*(x3 - x6)*(y2 - y3))/JE2r4);
    a(1,5)=-(((lambda + 2*mu)*(x2 - x3)*(x3 - x4) + 2*(lambda + mu)*(y2 - y3)*(y3 - y4))/JE1r3);
    a(1,6)=(lambda*(x3 - x4)*(y2 - y3))/JE1r3;
    a(2,1)=a(1,2);
    a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2)/JE1r3 + ...
           (2*(lambda + mu)*(x3 - x6)^2 + (lambda + 2*mu)*(y3 - y6)^2)/JE2r4;
    a(2,3)=-((lambda*(x2 - x3)*(y3 - y6))/JE2r4);
    a(2,4)=(2*(lambda + mu)*(x2 - x3)*(x3 - x6) + (lambda + 2*mu)*(y2 - y3)*(y3 - y6))/JE2r4;
    a(2,5)=(lambda*(x2 - x3)*(y3 - y4))/JE1r3;
    a(2,6)=(-2*(lambda + mu)*(x2 - x3)*(x3 - x4) - (lambda + 2*mu)*(y2 - y3)*(y3 - y4))/JE1r3;
    a(3,1)=a(1,3);
    a(3,2)=a(2,3);
    a(3,3)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/JE2r4;
    a(3,4)=-((lambda*(x2 - x3)*(y2 - y3))/JE2r4);
%     a(3,5)=0;
%     a(3,6)=0;
    a(4,1)= a(1,4);
    a(4,2)= a(2,4);
    a(4,3)= a(3,4);
    a(4,4)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/JE2r4;
%     a(4,5)=0;
%     a(4,6)=0;
    a(5,1)=a(1,5);
    a(5,2)=a(2,5);
%     a(5,3)=0;
%     a(5,4)=0;
    a(5,5)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2)/JE1r3;
    a(5,6)=-((lambda*(x2 - x3)*(y2 - y3))/JE1r3);
    a(6,1)=a(1,6);
    a(6,2)=a(2,6);
%     a(6,3)=0;
%     a(6,4)=0;
    a(6,5)=a(5,6);
    a(6,6)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2)/JE1r3;
    a=(1/(16.*mu*(lambda + mu)))*a;
    
    ind=ld+1:ld+vdim;
    matrix(ind,ind)=a;
%         tf = issymmetric(a)
%         d = eig(a);
%         isposdef = all(d > 0)
    ld=ld+vdim;
end

% North-East corner node (j=N+1)
i=N+1;

vdim=4;
a=zeros(vdim,vdim);

% x1=x(i - 1, j - 1);  
% y1=y(i - 1, j - 1);            
x2=x(i, j - 1);
y2=y(i, j - 1);
x3=x(i, j);
y3=y(i, j);
x4=x(i - 1, j);
y4=y(i - 1, j);

Jer3=abs(x4*(y2 - y3) + x2*(y3 - y4) + x3*(-y2 + y4));
denom=16*Jer3*mu*(lambda + mu);

a(1,1)=((lambda + 2*mu)*(x3 - x4)^2 + 2*(lambda + mu)*(y3 - y4)^2);
a(1,2)=-((lambda*(x3 - x4)*(y3 - y4)));
a(1,3)=-(((lambda + 2*mu)*(x2 - x3)*(x3 - x4) + 2*(lambda + mu)*(y2 - y3)*(y3 - y4)));
a(1,4)=(lambda*(x3 - x4)*(y2 - y3));
a(2,1)=a(1,2);
a(2,2)=(2*(lambda + mu)*(x3 - x4)^2 + (lambda + 2*mu)*(y3 - y4)^2);
a(2,3)=(lambda*(x2 - x3)*(y3 - y4));
a(2,4)=-((2*(lambda + mu)*(x2 - x3)*(x3 - x4) + (lambda + 2*mu)*(y2 - y3)*(y3 - y4)));
a(3,1)=a(1,3);
a(3,2)=a(2,3);
a(3,3)=((lambda + 2*mu)*(x2 - x3)^2 + 2*(lambda + mu)*(y2 - y3)^2);
a(3,4)=-((lambda*(x2 - x3)*(y2 - y3)));
a(4,1)=a(1,4);
a(4,2)=a(2,4);
a(4,3)=a(3,4);
a(4,4)=(2*(lambda + mu)*(x2 - x3)^2 + (lambda + 2*mu)*(y2 - y3)^2);
a=(1/denom)*a;

ind=ld+1:ld+vdim;
matrix(ind,ind)=a;
%     tf = issymmetric(a)
%     d = eig(a);
%     isposdef = all(d > 0)
end
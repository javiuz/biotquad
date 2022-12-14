function matrix=build_matrix_p_tau(alpha,lambda,mu)

% Creamos la matriz A_sigma_p^transp

global NN x y

N=NN;

matrix=sparse(8*N*(N+1),N*N);

coef=alpha/(8*(lambda+mu));

% South-West corner node
i=1;
j=1;

vdim=4;

a=zeros(vdim,1);

x1=x(i,j);  
y1=y(i,j);            
x2=x(i+1,j);
y2=y(i+1,j);
% x3=x(i+1,j+1);
% y3=y(i+1,j+1);
x4=x(i,j+1);
y4=y(i,j+1);

a(1,1)=(-x1 + x4);
a(2,1)=(-y1 + y4);
a(3,1)=(-x1 + x2);
a(4,1)=(-y1 + y2);

indr=1:vdim;
indc=i+(j-1)*N;
matrix(indr,indc)=a;
ld=vdim;

% South nodes (j=1)
vdim=6;
a=zeros(vdim,2);

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

%     a(1,1)=0;
%     a(2,1)=0;
    a(3,1)=(-x1 + x2);
    a(4,1)=(-y1 + y2);
    a(5,1)=(-x2 + x3);
    a(6,1)=(-y2 + y3);
    a(1,2)=(-x2 + x3);
    a(2,2)=(-y2 + y3);
    a(3,2)=(-x2 + x5);
    a(4,2)=(-y2 + y5);
%     a(5,2)=0;
%     a(6,2)=0;

    indr=ld+1:ld+vdim;
    indc=i+(j-1)*N;
    matrix(indr,indc-1:indc)=a;
    ld=ld+vdim;
end

% South-East corner node (j=1)
i=N+1;

vdim=4;
a=zeros(vdim,1);

x1=x(i-1,j);  
y1=y(i-1,j);            
x2=x(i,j);
y2=y(i,j);
x3=x(i,j+1);
y3=y(i,j+1);
% x4=x(i-1,j+1);
% y4=y(i-1,j+1);

a(1,1)=(-x1 + x2);
a(2,1)=(-y1 + y2);
a(3,1)=(-x2 + x3);
a(4,1)=(-y2 + y3);

indr=ld+1:ld+vdim;
indc=i-1+(j-1)*N;
matrix(indr,indc)=a;
ld=ld+vdim;


for j=2:N

    % West nodes
    i=1;
    
    vdim=6;
    a=zeros(vdim,2);  
    
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
    
    a(1,1)=(x3 - x4);
    a(2,1)=(y3 - y4);
    a(3,1)=(-x1 + x4); 
    a(4,1)=(-y1 + y4);
%     a(5,1)=0;
%     a(6,1)=0;
%     a(1,2)=0;
%     a(2,2)=0;
    a(3,2)=(-x4 + x6); 
    a(4,2)=(-y4 + y6);
    a(5,2)=(-x4 + x3);
    a(6,2)=(-y4 + y3);
    
    indr=ld+1:ld+vdim;
    indc=i+(j-1)*N;
    matrix(indr,[indc-N,indc])=a;
    ld=ld+vdim;
    
    % Central nodes 
    vdim=8;
    a=zeros(vdim,4);
    
    for i=2:N
        
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

    a(1,1)=(x3 - x4);
    a(2,1)=(y3 - y4);
%     a(3,1)=0;
%     a(4,1)=0;
%     a(5,1)=0;
%     a(6,1)=0;
    a(7,1)=(-x2 + x3);
    a(8,1)=(-y2 + y3);
    a(1,2)=(x6 - x3);
    a(2,2)=(y6 - y3);
    a(3,2)=(-x2 + x3);
    a(4,2)=(-y2 + y3);
%     a(5,2)=0;
%     a(6,2)=0;
%     a(7,2)=0;
%     a(8,2)=0;
%     a(1,3)=0;
%     a(2,3)=0;
    a(3,3)=(-x3 + x8);
    a(4,3)=(-y3 + y8);
    a(5,3)=(-x3 + x6);
    a(6,3)=(-y3 + y6);
%     a(7,3)=0;
%     a(8,3)=0;
%     a(1,4)=0;
%     a(2,4)=0;
%     a(3,4)=0;
%     a(4,4)=0;
    a(5,4)=(-x4 + x3);
    a(6,4)=(-y4 + y3);
    a(7,4)=(-x3 + x8);
    a(8,4)=(-y3 + y8);
    
    indr=ld+1:ld+vdim;
    indc=i+(j-1)*N;
    
    matrix(indr,indc-(N+1))=a(:,1);
    matrix(indr,indc-N)=a(:,2);
    matrix(indr,indc)=a(:,3);
    matrix(indr,indc-1)=a(:,4);
    
    ld=ld+vdim;
    end
    
    % East nodes
    i=N+1;
    
    vdim=6;
    a=zeros(vdim,2);
    
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

    a(1,1)=(x3 - x4);
    a(2,1)=(y3 - y4);
%     a(3,1)=0;
%     a(4,1)=0;
    a(5,1)=(-x2 + x3);
    a(6,1)=(-y2 + y3);
%     a(1,2)=0;
%     a(2,2)=0;
    a(3,2)=(-x4 + x3);
    a(4,2)=(-y4 + y3);
    a(5,2)=(-x3 + x5);
    a(6,2)=(-y3 + y5);
    
    indr=ld+1:ld+vdim;
    indc=i-1+(j-1)*N;
    matrix(indr,[indc-N,indc])=a;
    ld=ld+vdim;        
end

% North-West corner node
i=1;
j=N+1;

vdim=4;
a=zeros(vdim,1);

x1=x(i,j-1);
y1=y(i,j-1);           
% x2=x(i+1,j-1);
% y2=y(i+1,j-1);
x3=x(i+1,j);
y3=y(i+1,j);
x4=x(i,j);
y4=y(i,j);

a(1,1)=(x3 - x4);
a(2,1)=(y3 - y4);
a(3,1)=(-x1 + x4);
a(4,1)=(-y1 + y4);

indr=ld+1:ld+vdim;
indc=i+(j-2)*N;
matrix(indr,indc)=a;
ld=ld+vdim;

% North nodes (j=N+1)
vdim=6;
a=zeros(vdim,2);

for i=2:N
    
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
    
    a(1,1)=(x3 - x4);
    a(2,1)=(y3 - y4);
%     a(3,1)=0;
%     a(4,1)=0;
    a(5,1)=(-x2 + x3);
    a(6,1)=(-y2 + y3);
    a(1,2)=(x6 - x3);
    a(2,2)=(y6 - y3);
    a(3,2)=(-x2 + x3);
    a(4,2)=(-y2 + y3);
%     a(5,2)=0;
%     a(6,2)=0;

    indr=ld+1:ld+vdim;
    indc=i+(j-2)*N;
    matrix(indr,indc-1:indc)=a;
    ld=ld+vdim;
end

% North-East corner node (j=N+1)
i=N+1;

vdim=4;
a=zeros(vdim,1);

% x1=x(i-1,j-1);  
% y1=y(i-1,j-1);            
x2=x(i,j-1);
y2=y(i,j-1);
x3=x(i,j);
y3=y(i,j);
x4=x(i-1,j);
y4=y(i-1,j);

a(1,1)=(x3 - x4);
a(2,1)=(y3 - y4);
a(3,1)=(-x2 + x3);
a(4,1)=(-y2 + y3);

indr=ld+1:ld+vdim;
indc=i-1+(j-2)*N;
matrix(indr,indc)=a;

matrix=coef*matrix;
end
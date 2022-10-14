function AspT=build_matrix_sig_p

% Creamos la matriz A_sigma_p^transp

global NN x y alpha lambda mu

N=NN;

AspT=sparse(8*N*(N+1),N*N);

coef_d=alpha/(16*mu*(lambda + mu));

% South-West corner node
i=1;
j=1;
ind1p=(i+(j-1)*N);
vdim=4;

x1=x(i,j);  
y1=y(i,j);            
x2=x(i+1,j);
y2=y(i+1,j);
% x3=x(i+1,j+1);
% y3=y(i+1,j+1);
x4=x(i,j+1);
y4=y(i,j+1);

Jer1=abs(x4*(y1 - y2) + x1*(y2 - y4) + x2*(-y1 + y4));
denom_d=coef_d/Jer1;

% Matriz local (A sigma, p)^t
d=zeros(vdim,1);
d(1,1)=(lambda + 2*mu)*(x1 - x2)*(x1 - x4) - lambda*(x1 - x4)*(y1 - y4) + 2*(lambda + mu)*(y1 - y2)*(y1 - y4);
d(2,1)=2*(lambda + mu)*(x1 - x4)^2 + lambda*(-x1 + x2)*(y1 - y4) + (lambda + 2*mu)*(y1 - y4)^2;
d(3,1)=(lambda + 2*mu)*(x1 - x2)^2 + 2*(lambda + mu)*(y1 - y2)^2 + lambda*(-x1 + x2)*(y1 - y4);
d(4,1)=2*(lambda + mu)*(x1 - x2)*(x1 - x4) + lambda*(-x1 + x2)*(y1 - y2) + (lambda + 2*mu)*(y1 - y2)*(y1 - y4);
d=denom_d*d;

indr=1:vdim;
AspT(indr,ind1p)=d;
ld=vdim;

% South nodes (j=1)
vdim=6;

d=zeros(vdim,2);

for i=2:N
    ind1p=(i-1+(j-1)*N);
    ind2p=(i+(j-1)*N);
    
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
    denom_d1=coef_d/JE1r2;
    JE2r1=abs(x5*(-y2 + y3) + x3*(y2 - y5) + x2*(-y3 + y5));
    denom_d2=coef_d/JE2r1;

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

    indr=ld+1:ld+vdim;
    AspT(indr,ind1p:ind2p)=d;
    ld=ld+vdim;
end

% South-East corner node (j=1)
i=N+1;
ind1p=(i-1+(j-1)*N);
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
denom_d=coef_d/Jer2;

d=zeros(vdim,1);
d(1,1)=(lambda + 2*mu)*(x1 - x2)^2 + 2*(lambda + mu)*(y1 - y2)^2 + lambda*(-x1 + x2)*(y2 - y3);
d(2,1)=2*(lambda + mu)*(x1 - x2)*(x2 - x3) + lambda*(-x1 + x2)*(y1 - y2) + (lambda + 2*mu)*(y1 - y2)*(y2 - y3);
d(3,1)=(lambda + 2*mu)*(x1 - x2)*(x2 - x3) + lambda*(-x2 + x3)*(y2 - y3) + 2*(lambda + mu)*(y1 - y2)*(y2 - y3);
d(4,1)=2*(lambda + mu)*(x2 - x3)^2 - lambda*(x1 - x2)*(y2 - y3) + (lambda + 2*mu)*(y2 - y3)^2;
d=denom_d*d;

indr=ld+1:ld+vdim;
AspT(indr,ind1p)=d;
ld=ld+vdim;


for j=2:N

    % West nodes
    i=1;
    
    ind1p=(i+(j-2)*N);
    ind2p=(i+(j-1)*N);
    
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
    denom_d1=coef_d/JE1r4;
    JE2r1=abs(x6*(-y3 + y4) + x4*(y3 - y6) + x3*(-y4 + y6));
    denom_d2=coef_d/JE2r1;
    
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
    
    indr=ld+1:ld+vdim;
    AspT(indr,[ind1p,ind2p])=d;
    ld=ld+vdim;
    
    % Central nodes 
    vdim=8;
    
    d=zeros(vdim,4);
    
    for i=2:N
        ind1p=(i-1+(j-2)*N);
        ind2p=(i+(j-2)*N);
        ind3p=(i-1+(j-1)*N);
        ind4p=(i+(j-1)*N);
        
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
    denom_d1=coef_d/JE1r3;
    JE2r4=abs(x6*(-y2 + y3) + x3*(y2 - y6) + x2*(-y3 + y6));
    denom_d2=coef_d/JE2r4;
    JE3r1=abs(x8*(y3 - y6) + x3*(y6 - y8) + x6*(-y3 + y8));
    denom_d3=coef_d/JE3r1;
    JE4r2=abs(x8*(-y3 + y4) + x4*(y3 - y8) + x3*(-y4 + y8));
    denom_d4=coef_d/JE4r2;
    
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
    
    indr=ld+1:ld+vdim;
    AspT(indr,[ind1p:ind2p,ind3p:ind4p])=d;   
    ld=ld+vdim;
    end
    
    % East nodes
    i=N+1;
    ind1p=(i-1+(j-2)*N);
    ind2p=(i-1+(j-1)*N);
    
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
    denom_d1=coef_d/JE1r3;
    JE2r2=abs(x5*(-y3 + y4) + x4*(y3 - y5) + x3*(-y4 + y5));
    denom_d2=coef_d/JE2r2;
    
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
    
    indr=ld+1:ld+vdim;
    AspT(indr,[ind1p,ind2p])=d;
    ld=ld+vdim;        
end

% North-West corner node
i=1;
j=N+1;
ind1p=(i+(j-2)*N);
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
denom_d=coef_d/Jer4;

d=zeros(vdim,1);
d(1,1)=(lambda + 2*mu)*(x3 - x4)^2 + lambda*(x3 - x4)*(y1 - y4) + 2*(lambda + mu)*(y3 - y4)^2;
d(2,1)=-2*(lambda + mu)*(x1 - x4)*(x3 - x4) - lambda*(x3 - x4)*(y3 - y4) - (lambda + 2*mu)*(y1 - y4)*(y3 - y4);
d(3,1)=-((lambda + 2*mu)*(x1 - x4)*(x3 - x4)) - lambda*(x1 - x4)*(y1 - y4) - 2*(lambda + mu)*(y1 - y4)*(y3 - y4);
d(4,1)=2*(lambda + mu)*(x1 - x4)^2 + lambda*(x3 - x4)*(y1 - y4) + (lambda + 2*mu)*(y1 - y4)^2;
d=denom_d*d;

indr=ld+1:ld+vdim;
AspT(indr,ind1p)=d;
ld=ld+vdim;

% North nodes (j=N+1)
vdim=6;

d=zeros(vdim,2);

for i=2:N
    ind1p=(i-1+(j-2)*N);
    ind2p=(i+(j-2)*N);
    
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
    denom_d1=coef_d/JE1r3;
    JE2r4=abs(x6*(-y2 + y3) + x3*(y2 - y6) + x2*(-y3 + y6));
    denom_d2=coef_d/JE2r4;
    
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

    indr=ld+1:ld+vdim;
    AspT(indr,ind1p:ind2p)=d;
    ld=ld+vdim;
end

% North-East corner node (j=N+1)
i=N+1;
ind1p=(i-1+(j-2)*N);
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
denom_d=coef_d/Jer3;

d=zeros(vdim,1);
d(1,1)=(lambda + 2*mu)*(x3 - x4)^2 + lambda*(x3 - x4)*(y2 - y3) + 2*(lambda + mu)*(y3 - y4)^2;
d(2,1)=-2*(lambda + mu)*(x2 - x3)*(x3 - x4) + lambda*(-x3 + x4)*(y3 - y4) + (lambda + 2*mu)*(-y2 + y3)*(y3 - y4);
d(3,1)=-((lambda + 2*mu)*(x2 - x3)*(x3 - x4)) - lambda*(x2 - x3)*(y2 - y3) - 2*(lambda + mu)*(y2 - y3)*(y3 - y4);
d(4,1)=2*(lambda + mu)*(x2 - x3)^2 + lambda*(x3 - x4)*(y2 - y3) + (lambda + 2*mu)*(y2 - y3)^2;
d=denom_d*d;

indr=ld+1:ld+vdim;
AspT(indr,ind1p)=d;

end
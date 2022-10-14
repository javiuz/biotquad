function matrix=build_matrix_sig_gamma

% Creamos la matriz A_sigma_gamma^transp

global NN x y

N=NN;

matrix=sparse(8*N*(N+1),(N+1)*(N+1));
cont=1;

% South-West corner node
i=1;
j=1;

vdim=4;

x1=x(i,j);  
y1=y(i,j);            
x2=x(i+1,j);
y2=y(i+1,j);
% x3=x(i+1,j+1);
% y3=y(i+1,j+1);
x4=x(i,j+1);
y4=y(i,j+1);

% Matriz local (A sigma, gamma)^t
c=[y4-y1 x1-x4 y2-y1 x1-x2]';

ind=1:vdim;

matrix(ind,cont)=c;
ld=vdim;
cont=cont+1;

% South nodes (j=1)
vdim=6;

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
    
    c=[y3-y2 x2-x3 y5-y1 x1-x5 y3-y2 x2-x3]';
    
    ind=ld+1:ld+vdim;
    matrix(ind,cont)=c;

    ld=ld+vdim;
    cont=cont+1;
end

% South-East corner node (j=1)
i=N+1;

vdim=4;

x1=x(i-1,j);  
y1=y(i-1,j);            
x2=x(i,j);
y2=y(i,j);
x3=x(i,j+1);
y3=y(i,j+1);
% x4=x(i-1,j+1);
% y4=y(i-1,j+1);

c=[y2-y1 x1-x2 y3-y2 x2-x3]';

ind=ld+1:ld+vdim;
matrix(ind,cont)=c;
ld=ld+vdim;
cont=cont+1;

for j=2:N
    % West nodes
    i=1;
    
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
    
    c=[y3-y4 x4-x3 y6-y1 x1-x6 y3-y4 x4-x3]';
    
    ind=ld+1:ld+vdim;
    matrix(ind,cont)=c;
    ld=ld+vdim;
    cont=cont+1;
    
     
     % Central nodes 
    vdim=8;
    
    for i=2:N
%     x1=x(i-1,j-1);
%     y1=y(i-1,j-1);
    x2=x(i,j-1);
    y2=y(i,j-1);
%     x3=x(i,j);
%     y3=y(i,j);
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

    c=[y6-y4 x4-x6 y8-y2 x2-x8 y6-y4 x4-x6 y8-y2 x2-x8]';
    
    ind=ld+1:ld+vdim;
    matrix(ind,cont)=c;
    ld=ld+vdim;
    cont=cont+1;
    end
    
    % East nodes
    i=N+1;
  
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

    c=[y3-y4 x4-x3 y3-y4 x4-x3 y5-y2 x2-x5]';
    
    ind=ld+1:ld+vdim;
    matrix(ind,cont)=c;
    ld=ld+vdim;
    cont=cont+1;
end

% North-West corner node
i=1;
j=N+1;

vdim=4;

x1=x(i,j-1);
y1=y(i,j-1);           
% x2=x(i+1,j-1);
% y2=y(i+1,j-1);
x3=x(i+1,j);
y3=y(i+1,j);
x4=x(i,j);
y4=y(i,j);

c=[y3-y4 x4-x3 y4-y1 x1-x4]';

ind=ld+1:ld+vdim;
matrix(ind,cont)=c;
ld=ld+vdim;
cont=cont+1;

% North nodes (j=N+1)
vdim=6;

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
    
    c=[y6-y4 x4-x6 y3-y2 x2-x3 y3-y2 x2-x3]';
    
    ind=ld+1:ld+vdim;
    matrix(ind,cont)=c;

    ld=ld+vdim;
    cont=cont+1;
end

% North-East corner node (j=N+1)
i=N+1;

vdim=4;

% x1=x(i-1,j-1);  
% y1=y(i-1,j-1);            
x2=x(i,j-1);
y2=y(i,j-1);
x3=x(i,j);
y3=y(i,j);
x4=x(i-1,j);
y4=y(i-1,j);

c=[y3-y4 x4-x3 y3-y2 x2-x3]';

ind=ld+1:ld+vdim;
matrix(ind,cont)=c;

matrix=1/4*matrix;
end
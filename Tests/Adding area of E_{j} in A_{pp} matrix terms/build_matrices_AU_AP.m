function [matrix_u,matrix_p]=build_matrices_AU_AP

global NN x y

N=NN;

matrix_u=sparse(8*N*(N+1),2*N*N);
matrix_p=sparse(4*N*(N+1),N*N);

% South-West corner node
i=1;
j=1;

vdim=4;
vdim2=vdim/2;
a=zeros(vdim,2);
b=zeros(vdim2,1);

a(1,1)=1/(sqrt((x(i+1,j)-x(i,j))^2+(y(i+1,j)-y(i,j))^2));
a(2,2)=a(1,1);
a(3,1)=1/(sqrt((x(i,j+1)-x(i,j))^2+(y(i,j+1)-y(i,j))^2));
a(4,2)=a(3,1);

b(1,1)=a(1,1);
b(2,1)=a(3,1);

ind2u=(i+(j-1)*N)*2;
ind1u=ind2u-1;
ind1p=ind2u/2;

indru=1:vdim;
indrp=1:vdim2;

matrix_u(indru,ind1u:ind2u)=a;
matrix_p(indrp,ind1p)=b;
ld1=vdim;
ld2=vdim2;

% South nodes (j=1)
vdim=6;
vdim2=vdim/2;
a1=zeros(vdim,2);
a2=a1;
b1=zeros(vdim2,1);
b2=b1;

for i=2:N
    ind2u=(i-1+(j-1)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i+(j-1)*N)*2;
    ind3u=ind4u-1;  
    
    ind1p=ind2u/2;
    ind2p=ind4u/2;
    
    a1(3,1)=1/(sqrt((x(i,j+1)-x(i,j))^2+(y(i,j+1)-y(i,j))^2));
    a1(4,2)=a1(3,1);
    a1(5,1)=1/(sqrt((x(i-1,j)-x(i,j))^2+(y(i-1,j)-y(i,j))^2));
    a1(6,2)=a1(5,1);
    
    b1(2,1)=a1(3,1);
    b1(3,1)=a1(5,1);

    a2(1,1)=1/(sqrt((x(i+1,j)-x(i,j))^2+(y(i+1,j)-y(i,j))^2));
    a2(2,2)=a2(1,1);
    a2(3,1)=a1(3,1);
    a2(4,2)=a1(4,2);
    
    b2(1,1)=a2(1,1);
    b2(2,1)=b1(2,1);
    
    indru=ld1+1:ld1+vdim;
    indrp=ld2+1:ld2+vdim2;
    
    matrix_u(indru,ind1u:ind2u)=a1;
    matrix_u(indru,ind3u:ind4u)=a2;
    matrix_p(indrp,ind1p)=b1;
    matrix_p(indrp,ind2p)=b2;
    ld1=ld1+vdim;
    ld2=ld2+vdim2;
end

% South-East corner node (j=1)
i=N+1;

vdim=4;
vdim2=vdim/2;
a=zeros(vdim,2);
b=zeros(vdim2,1);

a(1,1)=1/(sqrt((x(i,j+1)-x(i,j))^2+(y(i,j+1)-y(i,j))^2));
a(2,2)=a(1,1);
a(3,1)=1/(sqrt((x(i-1,j)-x(i,j))^2+(y(i-1,j)-y(i,j))^2));
a(4,2)=a(3,1);

b(1,1)=a(1,1);
b(2,1)=a(3,1);

ind2u=(i-1+(j-1)*N)*2;
ind1u=ind2u-1;
ind1p=ind2u/2;

indru=ld1+1:ld1+vdim;
indrp=ld2+1:ld2+vdim2;
matrix_u(indru,ind1u:ind2u)=a;
matrix_p(indrp,ind1p)=b;
ld1=ld1+vdim;
ld2=ld2+vdim2;

for j=2:N
    % West nodes
    i=1;
    
    vdim=6;
    vdim2=vdim/2;
    a1=zeros(vdim,2);
    a2=a1; 
    b1=zeros(vdim2,1);
    b2=b1;
    
    a1(1,1)=1/(sqrt((x(i,j-1)-x(i,j))^2+(y(i,j-1)-y(i,j))^2));
    a1(2,2)=a1(1,1);
    a1(3,1)=1/(sqrt((x(i+1,j)-x(i,j))^2+(y(i+1,j)-y(i,j))^2));
    a1(4,2)=a1(3,1);
    
    b1(1,1)=a1(1,1);
    b1(2,1)=a1(3,1);
    
    a2(3,1)=a1(3,1);
    a2(4,2)=a1(4,2);
    a2(5,1)=1/(sqrt((x(i,j+1)-x(i,j))^2+(y(i,j+1)-y(i,j))^2));
    a2(6,2)=a2(5,1);
    
    b2(2,1)=b1(2,1);
    b2(3,1)=a2(5,1);
    
    ind2u=(i+(j-2)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i+(j-1)*N)*2;
    ind3u=ind4u-1;
    
    ind1p=ind2u/2;
    ind2p=ind4u/2;
    
    indru=ld1+1:ld1+vdim;
    indrp=ld2+1:ld2+vdim2;
    matrix_u(indru,ind1u:ind2u)=a1;
    matrix_u(indru,ind3u:ind4u)=a2;
    matrix_p(indrp,ind1p)=b1;
    matrix_p(indrp,ind2p)=b2;
    ld1=ld1+vdim;
    ld2=ld2+vdim2;
    
    % Central nodes 
    vdim=8;
    vdim2=vdim/2;
    a1=zeros(vdim,2);
    a2=a1;
    a3=a1;
    a4=a1;
    b1=zeros(vdim2,1);
    b2=b1;
    b3=b1;
    b4=b1;
    
    for i=2:N
        ind2u=(i-1+(j-2)*N)*2;
        ind1u=ind2u-1;
        ind4u=(i+(j-2)*N)*2;
        ind3u=ind4u-1;
        ind6u=(i-1+(j-1)*N)*2;
        ind5u=ind6u-1;
        ind8u=(i+(j-1)*N)*2;
        ind7u=ind8u-1;
        
        ind1p=ind2u/2;
        ind2p=ind4u/2;
        ind3p=ind6u/2;
        ind4p=ind8u/2;
        
        a1(1,1)=1/(sqrt((x(i,j-1)-x(i,j))^2+(y(i,j-1)-y(i,j))^2));
        a1(2,2)=a1(1,1);
        a1(7,1)=1/(sqrt((x(i-1,j)-x(i,j))^2+(y(i-1,j)-y(i,j))^2));
        a1(8,2)=a1(7,1);
        
        b1(1,1)=a1(1,1);
        b1(4,1)=a1(7,1);
    
        a2(1,1)=a1(1,1);
        a2(2,2)=a1(2,2);
        a2(3,1)=1/(sqrt((x(i+1,j)-x(i,j))^2+(y(i+1,j)-y(i,j))^2));
        a2(4,2)=a2(3,1);
        
        b2(1,1)=b1(1,1);
        b2(2,1)=a2(3,1);
    
        a3(5,1)=1/(sqrt((x(i,j+1)-x(i,j))^2+(y(i,j+1)-y(i,j))^2));
        a3(6,2)=a3(5,1);
        a3(7,1)=a1(7,1);
        a3(8,2)=a1(8,2);
        
        b3(3,1)=a3(5,1);
        b3(4,1)=b1(4,1);
    
        a4(3,1)=a2(3,1);
        a4(4,2)=a2(4,2);
        a4(5,1)=a3(5,1);
        a4(6,2)=a3(6,2);
        
        b4(2,1)=b2(2,1);
        b4(3,1)=b3(3,1);
        
        indru=ld1+1:ld1+vdim;
        indrp=ld2+1:ld2+vdim2;
        matrix_u(indru,ind1u:ind2u)=a1;
        matrix_u(indru,ind3u:ind4u)=a2;
        matrix_u(indru,ind5u:ind6u)=a3;
        matrix_u(indru,ind7u:ind8u)=a4;
        matrix_p(indrp,ind1p)=b1;
        matrix_p(indrp,ind2p)=b2;
        matrix_p(indrp,ind3p)=b3;
        matrix_p(indrp,ind4p)=b4;
        ld1=ld1+vdim;
        ld2=ld2+vdim2;
    end
    
    % East nodes
    i=N+1;
  
    vdim=6;
    vdim2=vdim/2;
    a1=zeros(vdim,2);
    a2=a1;
    b1=zeros(vdim2,1);
    b2=b1;
    
    a1(1,1)=1/(sqrt((x(i,j-1)-x(i,j))^2+(y(i,j-1)-y(i,j))^2));
    a1(2,2)=a1(1,1);
    a1(5,1)=1/(sqrt((x(i-1,j)-x(i,j))^2+(y(i-1,j)-y(i,j))^2));
    a1(6,2)=a1(5,1);
    
    b1(1,1)=a1(1,1);
    b1(3,1)=a1(5,1);
    
    a2(3,1)=1/(sqrt((x(i,j+1)-x(i,j))^2+(y(i,j+1)-y(i,j))^2));
    a2(4,2)=a2(3,1);
    a2(5,1)=a1(5,1);
    a2(6,2)=a1(6,2);
    
    b2(2,1)=a2(3,1);
    b2(3,1)=b1(3,1);
    
    ind2u=(i-1+(j-2)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i-1+(j-1)*N)*2;
    ind3u=ind4u-1;
    
    ind1p=ind2u/2;
    ind2p=ind4u/2;
    
    indru=ld1+1:ld1+vdim;
    indrp=ld2+1:ld2+vdim2;
    matrix_u(indru,ind1u:ind2u)=a1;
    matrix_u(indru,ind3u:ind4u)=a2;
    matrix_p(indrp,ind1p)=b1;
    matrix_p(indrp,ind2p)=b2;
    ld1=ld1+vdim;
    ld2=ld2+vdim2;
end

% North-West corner node
i=1;
j=N+1;

vdim=4;
vdim2=vdim/2;
a=zeros(vdim,2);
b=zeros(vdim2,1);

a(1,1)=1/(sqrt((x(i,j-1)-x(i,j))^2+(y(i,j-1)-y(i,j))^2));
a(2,2)=a(1,1);
a(3,1)=1/(sqrt((x(i+1,j)-x(i,j))^2+(y(i+1,j)-y(i,j))^2));
a(4,2)=a(3,1);

b(1,1)=a(1,1);
b(2,1)=a(3,1);

ind2u=(i+(j-2)*N)*2;
ind1u=ind2u-1;
ind1p=ind2u/2;

indru=ld1+1:ld1+vdim;
indrp=ld2+1:ld2+vdim2;
matrix_u(indru,ind1u:ind2u)=a;
matrix_p(indrp,ind1p)=b;
ld1=ld1+vdim;
ld2=ld2+vdim2;

% North nodes
vdim=6;
vdim2=vdim/2;
a1=zeros(vdim,2);
a2=a1;
b1=zeros(vdim2,1);
b2=b1;

for i=2:N
    ind2u=(i-1+(j-2)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i+(j-2)*N)*2;
    ind3u=ind4u-1;
    
    ind1p=ind2u/2;
    ind2p=ind4u/2;
    
    a1(1,1)=1/(sqrt((x(i,j-1)-x(i,j))^2+(y(i,j-1)-y(i,j))^2));
    a1(2,2)=a1(1,1);
    a1(5,1)=1/(sqrt((x(i-1,j)-x(i,j))^2+(y(i-1,j)-y(i,j))^2));
    a1(6,2)=a1(5,1);
    
    b1(1,1)=a1(1,1);
    b1(3,1)=a1(5,1);

    a2(1,1)=a1(1,1);
    a2(2,2)=a1(2,2);
    a2(3,1)=1/(sqrt((x(i+1,j)-x(i,j))^2+(y(i+1,j)-y(i,j))^2));
    a2(4,2)=a2(3,1);
    
    b2(1,1)=b1(1,1);
    b2(2,1)=a2(3,1);
    
    indru=ld1+1:ld1+vdim;
    indrp=ld2+1:ld2+vdim2;
    matrix_u(indru,ind1u:ind2u)=a1;
    matrix_u(indru,ind3u:ind4u)=a2;
    matrix_p(indrp,ind1p)=b1;
    matrix_p(indrp,ind2p)=b2;
    ld1=ld1+vdim;
    ld2=ld2+vdim2;
end

% North-East corner node (j=N+1)
i=N+1;

vdim=4;
vdim2=vdim/2;
a=zeros(vdim,2);
b=zeros(vdim2,1);

a(1,1)=1/(sqrt((x(i,j-1)-x(i,j))^2+(y(i,j-1)-y(i,j))^2));
a(2,2)=a(1,1);
a(3,1)=1/(sqrt((x(i-1,j)-x(i,j))^2+(y(i-1,j)-y(i,j))^2));
a(4,2)=a(3,1);

b(1,1)=a(1,1);
b(2,1)=a(3,1);

ind2u=(i-1+(j-2)*N)*2;
ind1u=ind2u-1;
ind1p=ind2u/2;

indru=ld1+1:ld1+vdim;
indrp=ld2+1:ld2+vdim2;
matrix_u(indru,ind1u:ind2u)=a;
matrix_p(indrp,ind1p)=b;

matrix_u=matrix_u';
matrix_p=matrix_p';
end
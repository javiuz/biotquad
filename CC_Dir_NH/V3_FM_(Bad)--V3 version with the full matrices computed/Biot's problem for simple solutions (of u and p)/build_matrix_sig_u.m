function matrix=build_matrix_sig_u

% Creamos la matriz A_sigma_u^transp

global NN 

N=NN;

matrix=sparse(8*N*(N+1),2*N*N);

% South-West corner node
i=1;
j=1;
ind2u=(i+(j-1)*N)*2;
ind1u=ind2u-1;
vdim=4;

% Matriz local (A sigma, u)^t
b=[-1 0;0 -1;-1 0;0 -1];

indr=1:vdim;
matrix(indr,ind1u:ind2u)=b;
ld=vdim;

% South nodes (j=1)
vdim=6;

% Matrix b is the same for every South node
b=[0 0 -1 0;0 0 0 -1;1 0 -1 0;0 1 0 -1;-1 0 0 0;0 -1 0 0];

for i=2:N
    ind2u=(i-1+(j-1)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i+(j-1)*N)*2;
%     ind3u=ind4u-1;
    
    indr=ld+1:ld+vdim;
    matrix(indr,ind1u:ind4u)=b;
    ld=ld+vdim;
end

% South-East corner node (j=1)
i=N+1;
ind2u=(i-1+(j-1)*N)*2;
ind1u=ind2u-1;
vdim=4;

b=[1 0;0 1;-1 0;0 -1];

indr=ld+1:ld+vdim;
matrix(indr,ind1u:ind2u)=b;
ld=ld+vdim;

for j=2:N
    % West nodes
    i=1;
    ind2u=(i+(j-2)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i+(j-1)*N)*2;
    ind3u=ind4u-1;
    
    vdim=6;
    
    b1=zeros(vdim,2);
    b2=b1;
    
    b1(1,1)=-1;
    b1(2,2)=-1;
    b1(3,1)=1;
    b1(4,2)=1;
    
    b2(3,1)=-1;
    b2(4,2)=-1;
    b2(5,1)=-1;
    b2(6,2)=-1;
    
    indr=ld+1:ld+vdim;
    matrix(indr,ind1u:ind2u)=b1;
    matrix(indr,ind3u:ind4u)=b2;
    ld=ld+vdim;
    
    % Central nodes 
    vdim=8;
    
    b1=zeros(vdim,2);
    b2=b1;
    b3=b1;
    b4=b1;
    
    b1(1,1)=1;
    b1(2,2)=1;
    b1(7,1)=1;
    b1(8,2)=1;
    
    b2(1,1)=-1;
    b2(2,2)=-1;
    b2(3,1)=1;
    b2(4,2)=1;
    
    b3(5,1)=1;
    b3(6,2)=1;
    b3(7,1)=-1;
    b3(8,2)=-1;
    
    b4(3,1)=-1;
    b4(4,2)=-1;
    b4(5,1)=-1;
    b4(6,2)=-1;
    
    for i=2:N
        ind2u=(i-1+(j-2)*N)*2;
        ind1u=ind2u-1;
        ind4u=(i+(j-2)*N)*2;
        ind3u=ind4u-1;
        ind6u=(i-1+(j-1)*N)*2;
        ind5u=ind6u-1;
        ind8u=(i+(j-1)*N)*2;
        ind7u=ind8u-1;
        
        indr=ld+1:ld+vdim;
        matrix(indr,ind1u:ind2u)=b1;
        matrix(indr,ind3u:ind4u)=b2;
        matrix(indr,ind5u:ind6u)=b3;
        matrix(indr,ind7u:ind8u)=b4;
        ld=ld+vdim;
    end
    
    % East nodes
    i=N+1;
    ind2u=(i-1+(j-2)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i-1+(j-1)*N)*2;
    ind3u=ind4u-1;
    
    vdim=6;
    
    b1=zeros(vdim,2);
    b2=b1;
    
    b1(1,1)=1;
    b1(2,2)=1;
    b1(5,1)=1;
    b1(6,2)=1;
    
    b2(3,1)=1;
    b2(4,2)=1;
    b2(5,1)=-1;
    b2(6,2)=-1;
    
    indr=ld+1:ld+vdim;
    matrix(indr,ind1u:ind2u)=b1;
    matrix(indr,ind3u:ind4u)=b2;
    ld=ld+vdim;
end

% North-West corner node
i=1;
j=N+1;
ind2u=(i+(j-2)*N)*2;
ind1u=ind2u-1;
vdim=4;

b=[-1 0;0 -1;1 0;0 1];

indr=ld+1:ld+vdim;
matrix(indr,ind1u:ind2u)=b;
ld=ld+vdim;

% North nodes
vdim=6;

% Matrix b is the same for every North node
b=[1 0 -1 0;0 1 0 -1;0 0 1 0;0 0 0 1;1 0 0 0;0 1 0 0];

for i=2:N
    ind2u=(i-1+(j-2)*N)*2;
    ind1u=ind2u-1;
    ind4u=(i+(j-2)*N)*2;
%     ind3u=ind4u-1;
    
    indr=ld+1:ld+vdim;
    matrix(indr,ind1u:ind4u)=b;
    ld=ld+vdim;
end

% North-East corner node (j=N+1)
i=N+1;
ind2u=(i-1+(j-2)*N)*2;
ind1u=ind2u-1;
vdim=4;

b=[1 0;0 1;1 0;0 1];

ind2u=(i-1+(j-2)*N)*2;
ind1u=ind2u-1;

indr=ld+1:ld+vdim;
matrix(indr,ind1u:ind2u)=b;

matrix=1/2*matrix;
end
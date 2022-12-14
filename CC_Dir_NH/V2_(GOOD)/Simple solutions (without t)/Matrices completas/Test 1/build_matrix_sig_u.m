function matrix=build_matrix_sig_u

% Creamos la matriz A_sigma_u^transp

global NN 

N=NN;

matrix=sparse(8*N*(N+1),2*N*N);

% South-West corner node
i=1;
j=1;

vdim=4;
a=zeros(vdim,2);

a(1,1)=-1;
a(2,2)=-1;
a(3,1)=-1;
a(4,2)=-1;

ind2=(i+(j-1)*N)*2;
ind1=ind2-1;

indr=1:vdim;
matrix(indr,ind1:ind2)=a;
ld=vdim;

% South nodes (j=1)
vdim=6;
a1=zeros(vdim,2);
a2=a1;

a1(3,1)=1;
a1(4,2)=1;
a1(5,1)=-1;
a1(6,2)=-1;

a2(1,1)=-1;
a2(2,2)=-1;
a2(3,1)=-1;
a2(4,2)=-1;

for i=2:N
    ind2=(i-1+(j-1)*N)*2;
    ind1=ind2-1;
    ind4=(i+(j-1)*N)*2;
    ind3=ind4-1;
    
    indr=ld+1:ld+vdim;
    matrix(indr,ind1:ind2)=a1;
    matrix(indr,ind3:ind4)=a2;
    ld=ld+vdim;
end

% South-East corner node (j=1)
i=N+1;

vdim=4;
a=zeros(vdim,2);

a(1,1)=1;
a(2,2)=1;
a(3,1)=-1;
a(4,2)=-1;

ind2=(i-1+(j-1)*N)*2;
ind1=ind2-1;

indr=ld+1:ld+vdim;
matrix(indr,ind1:ind2)=a;
ld=ld+vdim;

for j=2:N
    % West nodes
    i=1;
    
    vdim=6;
    a1=zeros(vdim,2);
    a2=a1;
    
    a1(1,1)=-1;
    a1(2,2)=-1;
    a1(3,1)=1;
    a1(4,2)=1;
    
    a2(3,1)=-1;
    a2(4,2)=-1;
    a2(5,1)=-1;
    a2(6,2)=-1;
    
    ind2=(i+(j-2)*N)*2;
    ind1=ind2-1;
    ind4=(i+(j-1)*N)*2;
    ind3=ind4-1;
    
    indr=ld+1:ld+vdim;
    matrix(indr,ind1:ind2)=a1;
    matrix(indr,ind3:ind4)=a2;
    ld=ld+vdim;
    
    % Central nodes 
    vdim=8;
    a1=zeros(vdim,2);
    a2=a1;
    a3=a1;
    a4=a1;
    
    a1(1,1)=1;
    a1(2,2)=1;
    a1(7,1)=1;
    a1(8,2)=1;
    
    a2(1,1)=-1;
    a2(2,2)=-1;
    a2(3,1)=1;
    a2(4,2)=1;
    
    a3(5,1)=1;
    a3(6,2)=1;
    a3(7,1)=-1;
    a3(8,2)=-1;
    
    a4(3,1)=-1;
    a4(4,2)=-1;
    a4(5,1)=-1;
    a4(6,2)=-1;
    
    for i=2:N
        ind2=(i-1+(j-2)*N)*2;
        ind1=ind2-1;
        ind4=(i+(j-2)*N)*2;
        ind3=ind4-1;
        ind6=(i-1+(j-1)*N)*2;
        ind5=ind6-1;
        ind8=(i+(j-1)*N)*2;
        ind7=ind8-1;
        
        indr=ld+1:ld+vdim;
        matrix(indr,ind1:ind2)=a1;
        matrix(indr,ind3:ind4)=a2;
        matrix(indr,ind5:ind6)=a3;
        matrix(indr,ind7:ind8)=a4;
        ld=ld+vdim;
    end
    
    % East nodes
    i=N+1;
  
    vdim=6;
    a1=zeros(vdim,2);
    a2=a1;
    
    a1(1,1)=1;
    a1(2,2)=1;
    a1(5,1)=1;
    a1(6,2)=1;
    
    a2(3,1)=1;
    a2(4,2)=1;
    a2(5,1)=-1;
    a2(6,2)=-1;
    
    ind2=(i-1+(j-2)*N)*2;
    ind1=ind2-1;
    ind4=(i-1+(j-1)*N)*2;
    ind3=ind4-1;
    
    indr=ld+1:ld+vdim;
    matrix(indr,ind1:ind2)=a1;
    matrix(indr,ind3:ind4)=a2;
    ld=ld+vdim;
end

% North-West corner node
i=1;
j=N+1;

vdim=4;
a=zeros(vdim,2);

a(1,1)=-1;
a(2,2)=-1;
a(3,1)=1;
a(4,2)=1;

ind2=(i+(j-2)*N)*2;
ind1=ind2-1;

indr=ld+1:ld+vdim;
matrix(indr,ind1:ind2)=a;
ld=ld+vdim;

% North nodes
vdim=6;
a1=zeros(vdim,2);
a2=a1;

a1(1,1)=1;
a1(2,2)=1;
a1(5,1)=1;
a1(6,2)=1;

a2(1,1)=-1;
a2(2,2)=-1;
a2(3,1)=1;
a2(4,2)=1;

for i=2:N
    ind2=(i-1+(j-2)*N)*2;
    ind1=ind2-1;
    ind4=(i+(j-2)*N)*2;
    ind3=ind4-1;
    
    indr=ld+1:ld+vdim;
    matrix(indr,ind1:ind2)=a1;
    matrix(indr,ind3:ind4)=a2;
    ld=ld+vdim;
end

% North-East corner node (j=N+1)
i=N+1;

vdim=4;
a=zeros(vdim,2);

a(1,1)=1;
a(2,2)=1;
a(3,1)=1;
a(4,2)=1;

ind2=(i-1+(j-2)*N)*2;
ind1=ind2-1;

indr=ld+1:ld+vdim;
matrix(indr,ind1:ind2)=a;

matrix=1/2*matrix;
end
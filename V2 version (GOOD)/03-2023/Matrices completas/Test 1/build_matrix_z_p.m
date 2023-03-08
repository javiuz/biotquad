function matrix=build_matrix_z_p

% Creamos la matriz A_zp^transp

global NN

N=NN;

matrix=sparse(4*N*(N+1),N*N);

% South-West corner node
i=1;
j=1;
indc=i+(j-1)*N;

vdim=2;
indr=1:vdim;

b=[1;1];
matrix(indr,indc)=b;

ld=vdim;

% South nodes (j=1)
vdim=3;

% Matrix b is the same for every South node
b=[0 1;-1 1;1 0];

for i=2:N
    ind1=i-1+(j-1)*N;
    ind2=i+(j-1)*N;
    indc=ind1:ind2;
    indr=ld+1:ld+vdim;
    
    matrix(indr,indc)=b;
    ld=ld+vdim;
end

% South-East corner node (j=1)
i=N+1;
indc=i-1+(j-1)*N;

vdim=2;
indr=ld+1:ld+vdim;

b=[-1;1];

matrix(indr,indc)=b;
ld=ld+vdim;

for j=2:N

    % West nodes
    i=1;
    ind1=i+(j-2)*N;
    ind2=i+(j-1)*N;
    indc=[ind1,ind2];
    
    vdim=3;
    indr=ld+1:ld+vdim;
    
    b=[1 0;-1 1;0 1];
    
    matrix(indr,indc)=b;
    ld=ld+vdim;
    
    % Central nodes 
    vdim=4;
    
    % Matrix b is the same for every central node
    b=[-1 1 0 0;0 -1 0 1;0 0 -1 1;-1 0 1 0];
    
    for i=2:N
        
        ind1=i-1+(j-2)*N;
        ind2=i+(j-2)*N;
        ind3=i-1+(j-1)*N;
        ind4=i+(j-1)*N;
        
        indc=[ind1,ind2,ind3,ind4];
        indr=ld+1:ld+vdim;
        
        matrix(indr,indc)=b;
        ld=ld+vdim;
    end
    
    % East nodes
    i=N+1;
    ind1=i-1+(j-2)*N;
    ind2=i-1+(j-1)*N; 
    indc=[ind1,ind2];
    
    vdim=3;
    indr=ld+1:ld+vdim;
    
    b=[-1 0;0 -1;-1 1];
    
    matrix(indr,indc)=b;
    ld=ld+vdim;
end

% North-West corner node
i=1;
j=N+1;
indc=i+(j-2)*N;

vdim=2;
indr=ld+1:ld+vdim;

b=[1;-1];

matrix(indr,indc)=b;
ld=ld+vdim;

% North nodes (j=N+1)
vdim=3;

% Matrix b is the same for every North node
b=[-1 1;0 -1;-1 0];

for i=2:N
    ind1=i-1+(j-2)*N;
    ind2=i+(j-2)*N;
    indc=ind1:ind2;
    
    indr=ld+1:ld+vdim;
    matrix(indr,indc)=b;
    ld=ld+vdim;
end

% North-East corner node (j=N+1)
i=N+1;
indc=i-1+(j-2)*N;

vdim=2;
indr=ld+1:ld+vdim;

b=[-1;-1];

matrix(indr,indc)=b;

matrix=(1/2)*matrix;
end
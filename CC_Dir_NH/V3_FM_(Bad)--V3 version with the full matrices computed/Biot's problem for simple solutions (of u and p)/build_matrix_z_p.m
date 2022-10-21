function matrix=build_matrix_z_p

% Creamos la matriz A_zp^transp

global NN

N=NN;

matrix=sparse(4*N*(N+1),N*N);

% South-West corner node
i=1;
j=1;
indp=i+(j-1)*N;

vdim=2;
indr=1:vdim;

f=[1;1];
matrix(indr,indp)=f;

ld=vdim;

% South nodes (j=1)
vdim=3;

% Matrix f is the same for every South node
f=[0 1;-1 1;1 0];

for i=2:N
    ind1p=i-1+(j-1)*N;
    ind2p=i+(j-1)*N;
    indp=ind1p:ind2p;
    indr=ld+1:ld+vdim;
    
    matrix(indr,indp)=f;
    ld=ld+vdim;
end

% South-East corner node (j=1)
i=N+1;
indp=i-1+(j-1)*N;

vdim=2;
indr=ld+1:ld+vdim;

f=[-1;1];

matrix(indr,indp)=f;
ld=ld+vdim;

for j=2:N

    % West nodes
    i=1;
    ind1p=i+(j-2)*N;
    ind2p=i+(j-1)*N;
    indp=[ind1p,ind2p];
    
    vdim=3;
    indr=ld+1:ld+vdim;
    
    f=[1 0;-1 1;0 1];
    
    matrix(indr,indp)=f;
    ld=ld+vdim;
    
    % Central nodes 
    vdim=4;
    
    % Matrix f is the same for every central node
    f=[-1 1 0 0;0 -1 0 1;0 0 -1 1;-1 0 1 0];
    
    for i=2:N
        
        ind1p=i-1+(j-2)*N;
        ind2p=i+(j-2)*N;
        ind3p=i-1+(j-1)*N;
        ind4p=i+(j-1)*N;
        
        indp=[ind1p,ind2p,ind3p,ind4p];
        indr=ld+1:ld+vdim;
        
        matrix(indr,indp)=f;
        ld=ld+vdim;
    end
    
    % East nodes
    i=N+1;
    ind1p=i-1+(j-2)*N;
    ind2p=i-1+(j-1)*N; 
    indp=[ind1p,ind2p];
    
    vdim=3;
    indr=ld+1:ld+vdim;
    
    f=[-1 0;0 -1;-1 1];
    
    matrix(indr,indp)=f;
    ld=ld+vdim;
end

% North-West corner node
i=1;
j=N+1;
indp=i+(j-2)*N;

vdim=2;
indr=ld+1:ld+vdim;

f=[1;-1];

matrix(indr,indp)=f;
ld=ld+vdim;

% North nodes (j=N+1)
vdim=3;

% Matrix f is the same for every North node
f=[-1 1;0 -1;-1 0];

for i=2:N
    ind1p=i-1+(j-2)*N;
    ind2p=i+(j-2)*N;
    indp=ind1p:ind2p;
    
    indr=ld+1:ld+vdim;
    matrix(indr,indp)=f;
    ld=ld+vdim;
end

% North-East corner node (j=N+1)
i=N+1;
indp=i-1+(j-2)*N;

vdim=2;
indr=ld+1:ld+vdim;

f=[-1;-1];

matrix(indr,indp)=f;

matrix=(1/2)*matrix;
end
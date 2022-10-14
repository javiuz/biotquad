function matrix=build_matrix_z_z

% Creamos la matriz A_z_z

global NN

N=NN;

matrix=sparse(4*N*(N+1),4*N*(N+1));

% South-West corner node
i=1;
j=1;
vdim=2;

t=zeros(vdim,vdim);
t(1,1)=kinv(2,2,i,j,i,j);
t(1,2)=kinv(2,1,i,j,i,j);
t(2,1)=t(1,2);
t(2,2)=kinv(1,1,i,j,i,j);

ind=1:vdim;
matrix(ind,ind)=t;

ld=vdim;

% South nodes (j=1)
vdim=3;
t=zeros(vdim,vdim);

for i=2:N
    t(1,1)=kinv(2,2,i,j,i,j);
    t(1,2)=kinv(2,1,i,j,i,j);
    %a(1,3)=0
    t(2,1)=t(1,2);
    t(2,2)=kinv(1,1,i-1,j,i,j)+kinv(1,1,i,j,i,j);
    t(2,3)=kinv(1,2,i-1,j,i,j);
    %a(3,1)=0
    t(3,2)=t(2,3);
    t(3,3)=kinv(2,2,i-1,j,i,j);
    
    ind=ld+1:ld+vdim;
    matrix(ind,ind)=t;
    ld=ld+vdim;
end 

% South-East corner node (j=1)
i=N+1;
vdim=2;

t=zeros(vdim,vdim);

t(1,1)=kinv(1,1,i-1,j,i,j);
t(1,2)=kinv(1,2,i-1,j,i,j);
t(2,1)=t(1,2);
t(2,2)=kinv(2,2,i-1,j,i,j);

ind=ld+1:ld+vdim;
matrix(ind,ind)=t;
ld=ld+vdim;

for j=2:N

    % West nodes
    i=1;
    vdim=3;
    t=zeros(vdim,vdim);  
    
    t(1,1)=kinv(1,1,i,j-1,i,j);
    t(1,2)=kinv(1,2,i,j-1,i,j);
    %a(1,3)=0;
    t(2,1)=t(1,2);
    t(2,2)=kinv(2,2,i,j-1,i,j)+kinv(2,2,i,j,i,j);
    t(2,3)=kinv(2,1,i,j,i,j);
    %a(3,1)=0;
    t(3,2)=t(2,3);
    t(3,3)=kinv(1,1,i,j,i,j);
    
    ind=ld+1:ld+vdim;
    matrix(ind,ind)=t;
    ld=ld+vdim;
    
    % Central nodes 
    vdim=4;

    t=zeros(vdim,vdim);
    
    for i=2:N
        
        t(1,1)=kinv(1,1,i-1,j-1,i,j)+kinv(1,1,i,j-1,i,j);
        t(1,2)=kinv(1,2,i,j-1,i,j);
        %a(1,3)=0;
        t(1,4)=kinv(1,2,i-1,j-1,i,j);
        t(2,1)=t(1,2);
        t(2,2)=kinv(2,2,i,j-1,i,j)+kinv(2,2,i,j,i,j);
        t(2,3)=kinv(2,1,i,j,i,j);
        %a(2,4)=0;
        %a(3,1)=0;
        t(3,2)=t(2,3);
        t(3,3)=kinv(1,1,i,j,i,j)+kinv(1,1,i-1,j,i,j);
        t(3,4)=kinv(1,2,i-1,j,i,j);
        t(4,1)=t(1,4);
        %a(4,2)=0;
        t(4,3)=t(3,4);
        t(4,4)=kinv(2,2,i-1,j,i,j)+kinv(2,2,i-1,j-1,i,j);
        
        ind=ld+1:ld+vdim;
        matrix(ind,ind)=t;
        ld=ld+vdim;
    end
    
    % East nodes
    i=N+1;
    
    vdim=3;
    t=zeros(vdim,vdim);
    
    t(1,1)=kinv(1,1,i-1,j-1,i,j);
    %a(1,2)=0;
    t(1,3)=kinv(1,2,i-1,j-1,i,j);
    %a(2,1)=0;
    t(2,2)=kinv(1,1,i-1,j,i,j);
    t(2,3)=kinv(1,2,i-1,j,i,j);
    t(3,1)=t(1,3);
    t(3,2)=t(2,3);
    t(3,3)=kinv(2,2,i-1,j-1,i,j)+kinv(2,2,i-1,j,i,j);
    
    ind=ld+1:ld+vdim;
    matrix(ind,ind)=t;
    ld=ld+vdim;
end

% North-West corner node
i=1;
j=N+1;
vdim=2;
t=zeros(vdim,vdim);

t(1,1)=kinv(1,1,i,j-1,i,j);
t(1,2)=kinv(1,2,i,j-1,i,j);
t(2,1)=t(1,2);
t(2,2)=kinv(2,2,i,j-1,i,j);

ind=ld+1:ld+vdim;
matrix(ind,ind)=t;
ld=ld+vdim;

% North nodes (j=N+1)
vdim=3;
t=zeros(vdim,vdim);

for i=2:N
    
    t(1,1)=kinv(1,1,i-1,j-1,i,j)+kinv(1,1,i,j-1,i,j);
    t(1,2)=kinv(1,2,i,j-1,i,j);
    t(1,3)=kinv(1,2,i-1,j-1,i,j);
    t(2,1)=t(1,2);
    t(2,2)=kinv(2,2,i,j-1,i,j);
    %a(2,3)=0;
    t(3,1)=t(1,3);
    %a(3,2)=0;
    t(3,3)=kinv(2,2,i-1,j-1,i,j);
    
    ind=ld+1:ld+vdim;
    matrix(ind,ind)=t;
    ld=ld+vdim;
end

% North-East corner node (j=N+1)
i=N+1;

vdim=2;
t=zeros(vdim,vdim);

t(1,1)=kinv(1,1,i-1,j-1,i,j);
t(1,2)=kinv(1,2,i-1,j-1,i,j);
t(2,1)=t(1,2);
t(2,2)=kinv(2,2,i-1,j-1,i,j);

ind=ld+1:ld+vdim;
matrix(ind,ind)=t;

matrix=(1/4)*matrix;

end
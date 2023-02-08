function M3=build_flux_matrix_t0(delta_t)

global NN x y c0 

N=NN;

M3=sparse(N*N,N*N);

% South-West corner node
i=1;
j=1;
ind2u=(i+(j-1)*N)*2;
ind1p=ind2u/2;
vdim=4;
vdim2=vdim/2;

% Matriz local (A z, p)^T
f=(1/2)*[1;1];

% Matriz local (A z, z)^{-1}
t=zeros(vdim2,vdim2);
t(1,1)=kinv(2,2,i,j,i,j);
t(1,2)=kinv(2,1,i,j,i,j);
t(2,1)=t(1,2);
t(2,2)=kinv(1,1,i,j,i,j);
t=(1/4)*t;

% localm3=k +(delta_t*f')*(t\f);
localm3=(delta_t*f')*(t\f);

M3(ind1p,ind1p)=M3(ind1p,ind1p)+localm3;

% South nodes (j=1)
vdim=6;
vdim2=vdim/2;

% Matrix f is the same for every South node
f=(1/2)*[0 1;-1 1;1 0];

t=zeros(vdim2,vdim2);

for i=2:N
    ind2u=(i-1+(j-1)*N)*2;
    ind4u=(i+(j-1)*N)*2;
    
    ind1p=ind2u/2;
    ind2p=ind4u/2;

    t(1,1)=kinv(2,2,i,j,i,j);
    t(1,2)=kinv(2,1,i,j,i,j);
    %t(1,3)=0
    t(2,1)=t(1,2);
    t(2,2)=kinv(1,1,i-1,j,i,j)+kinv(1,1,i,j,i,j);
    t(2,3)=kinv(1,2,i-1,j,i,j);
    %t(3,1)=0
    t(3,2)=t(2,3);
    t(3,3)=kinv(2,2,i-1,j,i,j);
    t=(1/4)*t;
    
%     localm3=k +(delta_t*f')*(t\f);
    localm3=(delta_t*f')*(t\f);
    
    M3(ind1p:ind2p,ind1p:ind2p)=M3(ind1p:ind2p,ind1p:ind2p)+localm3;
    
end

% South-East corner node (j=1)
i=N+1;
ind2u=(i-1+(j-1)*N)*2;
ind1p=ind2u/2;
vdim=4;
vdim2=vdim/2;

f=(1/2)*[-1;1];

t=zeros(vdim2,vdim2);
t(1,1)=kinv(1,1,i-1,j,i,j);
t(1,2)=kinv(1,2,i-1,j,i,j);
t(2,1)=t(1,2);
t(2,2)=kinv(2,2,i-1,j,i,j);
t=(1/4)*t;

localm3=(delta_t*f')*(t\f);

M3(ind1p,ind1p)=M3(ind1p,ind1p)+localm3;

for j=2:N

    % West nodes
    i=1;
    ind2u=(i+(j-2)*N)*2;
    ind4u=(i+(j-1)*N)*2;
    
    ind1p=ind2u/2;
    ind2p=ind4u/2;
    
    vdim=6;
    vdim2=vdim/2;
    
    f=(1/2)*[1 0;-1 1;0 1];
    
    t=zeros(vdim2,vdim2);  
    t(1,1)=kinv(1,1,i,j-1,i,j);
    t(1,2)=kinv(1,2,i,j-1,i,j);
    %t(1,3)=0;
    t(2,1)=t(1,2);
    t(2,2)=kinv(2,2,i,j-1,i,j)+kinv(2,2,i,j,i,j);
    t(2,3)=kinv(2,1,i,j,i,j);
    %t(3,1)=0;
    t(3,2)=t(2,3);
    t(3,3)=kinv(1,1,i,j,i,j);
    t=(1/4)*t;
    
%     localm3=k +(delta_t*f')*(t\f);
    localm3=(delta_t*f')*(t\f);
   
    M3([ind1p,ind2p],[ind1p,ind2p])=M3([ind1p,ind2p],[ind1p,ind2p])+localm3;
    
    % Central nodes 
    vdim=8;
    vdim2=vdim/2;

    t=zeros(vdim2,vdim2);
    % Matrix f is the same for every central node
    f=(1/2)*[-1 1 0 0;0 -1 0 1;0 0 -1 1;-1 0 1 0];
         
    for i=2:N
        ind2u=(i-1+(j-2)*N)*2;
        ind4u=(i+(j-2)*N)*2;
        ind6u=(i-1+(j-1)*N)*2;
        ind8u=(i+(j-1)*N)*2;
        
        ind1p=ind2u/2;
        ind2p=ind4u/2;
        ind3p=ind6u/2;
        ind4p=ind8u/2;
        
    t(1,1)=kinv(1,1,i-1,j-1,i,j)+kinv(1,1,i,j-1,i,j);
    t(1,2)=kinv(1,2,i,j-1,i,j);
    %t(1,3)=0;
    t(1,4)=kinv(1,2,i-1,j-1,i,j);
    t(2,1)=t(1,2);
    t(2,2)=kinv(2,2,i,j-1,i,j)+kinv(2,2,i,j,i,j);
    t(2,3)=kinv(2,1,i,j,i,j);
    %t(2,4)=0;
    %t(3,1)=0;
    t(3,2)=t(2,3);
    t(3,3)=kinv(1,1,i,j,i,j)+kinv(1,1,i-1,j,i,j);
    t(3,4)=kinv(1,2,i-1,j,i,j);
    t(4,1)=t(1,4);
    %t(4,2)=0;
    t(4,3)=t(3,4);
    t(4,4)=kinv(2,2,i-1,j,i,j)+kinv(2,2,i-1,j-1,i,j);
    t=(1/4)*t;
    
%     localm3=k +(delta_t*f')*(t\f);
    localm3=(delta_t*f')*(t\f);
      
    M3([ind1p:ind2p,ind3p:ind4p],[ind1p:ind2p,ind3p:ind4p])=M3([ind1p:ind2p,ind3p:ind4p],[ind1p:ind2p,ind3p:ind4p])+localm3;
    
    end
    
    % East nodes
    i=N+1;
    ind2u=(i-1+(j-2)*N)*2;
    ind4u=(i-1+(j-1)*N)*2;

    ind1p=ind2u/2;
    ind2p=ind4u/2;
    
    vdim=6;
    vdim2=vdim/2;
    
    f=(1/2)*[-1 0;0 -1;-1 1];
    
    t=zeros(vdim2,vdim2);
    t(1,1)=kinv(1,1,i-1,j-1,i,j);
    %t(1,2)=0;
    t(1,3)=kinv(1,2,i-1,j-1,i,j);
    %t(2,1)=0;
    t(2,2)=kinv(1,1,i-1,j,i,j);
    t(2,3)=kinv(1,2,i-1,j,i,j);
    t(3,1)=t(1,3);
    t(3,2)=t(2,3);
    t(3,3)=kinv(2,2,i-1,j-1,i,j)+kinv(2,2,i-1,j,i,j);
    t=(1/4)*t;
    
%     k1=coef1_p;
%     k2=coef2_p;
%     k=[k1 0;0 k2];

    localm3=(delta_t*f')*(t\f);
    
    M3([ind1p,ind2p],[ind1p,ind2p])=M3([ind1p,ind2p],[ind1p,ind2p])+localm3;
end

% North-West corner node
i=1;
j=N+1;
ind2u=(i+(j-2)*N)*2;
ind1p=ind2u/2;
vdim=4;
vdim2=vdim/2;

f=(1/2)*[1;-1];

t=zeros(vdim2,vdim2);
t(1,1)=kinv(1,1,i,j-1,i,j);
t(1,2)=kinv(1,2,i,j-1,i,j);
t(2,1)=t(1,2);
t(2,2)=kinv(2,2,i,j-1,i,j);
t=(1/4)*t;

% k=coef_p;

localm3=(delta_t*f')*(t\f);

M3(ind1p,ind1p)=M3(ind1p,ind1p)+localm3;

% North nodes (j=N+1)
vdim=6;
vdim2=vdim/2;

% Matrix f is the same for every North node
f=(1/2)*[-1 1;0 -1;-1 0];

t=zeros(vdim2,vdim2);

for i=2:N
    ind2u=(i-1+(j-2)*N)*2;
    ind4u=(i+(j-2)*N)*2;
    
    ind1p=ind2u/2;
    ind2p=ind4u/2;
    
    t(1,1)=kinv(1,1,i-1,j-1,i,j)+kinv(1,1,i,j-1,i,j);
    t(1,2)=kinv(1,2,i,j-1,i,j);
    t(1,3)=kinv(1,2,i-1,j-1,i,j);
    t(2,1)=t(1,2);
    t(2,2)=kinv(2,2,i,j-1,i,j);
    %t(2,3)=0;
    t(3,1)=t(1,3);
    %t(3,2)=0;
    t(3,3)=kinv(2,2,i-1,j-1,i,j);
    t=(1/4)*t;

%     k1=coef1_p;
%     k2=coef2_p;
%     k=[k1 0;0 k2];

    localm3=(delta_t*f')*(t\f);

    M3(ind1p:ind2p,ind1p:ind2p)=M3(ind1p:ind2p,ind1p:ind2p)+localm3;
end

% North-East corner node (j=N+1)
i=N+1;
ind2u=(i-1+(j-2)*N)*2;
ind1p=ind2u/2;
vdim=4;
vdim2=vdim/2;

f=(1/2)*[-1;-1];

t=zeros(vdim2,vdim2);
t(1,1)=kinv(1,1,i-1,j-1,i,j);
t(1,2)=kinv(1,2,i-1,j-1,i,j);
t(2,1)=t(1,2);
t(2,2)=kinv(2,2,i-1,j-1,i,j);
t=(1/4)*t;

% k=coef_p;

localm3=(delta_t*f')*(t\f);

M3(ind1p,ind1p)=M3(ind1p,ind1p)+localm3;
return
end
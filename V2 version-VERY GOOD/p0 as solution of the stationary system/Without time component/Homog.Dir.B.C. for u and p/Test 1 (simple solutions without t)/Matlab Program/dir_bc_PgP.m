function gDp=dir_bc_PgP(delta_t,tt)

global NN 

N=NN;

gDp=zeros(N*N,1);

% South-West corner node
i=1;
j=1;
ind2u=(i+(j-1)*N)*2;
ind1p=ind2u/2;
vdim=4;
vdim2=vdim/2;

f=(1/2)*[1;1];

t=zeros(vdim2,vdim2);
t(1,1)=kinv(2,2,i,j,i,j);
t(1,2)=kinv(2,1,i,j,i,j);
t(2,1)=t(1,2);
t(2,2)=kinv(1,1,i,j,i,j);
t=(1/4)*t;

% CC.D en el nodo S-W, fórmula de cuadratura para el término Pg de p en:
    % frontera Sur
        Pgp_L1=int_dir_simpson(i,j,i+1,j,tt,3);
    % frontera Oeste
        Pgp_L4=int_dir_simpson(i,j,i,j+1,tt,3);
    
        Pgp=(1/2)*[Pgp_L1;Pgp_L4];  
        localv_p=((delta_t*f')*(t\Pgp));
        gDp(ind1p)=gDp(ind1p)+localv_p;

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
    
    % CC.D en los nodos Sur, fórmula de cuadratura para el término Pg de p en:
        % frontera Sur 1
            Pgp_L2=int_dir_simpson(i,j,i+1,j,tt,3); 
        % frontera Sur 2
            Pgp_L1=int_dir_simpson(i-1,j,i,j,tt,3); 
        Pgp=(1/2)*[Pgp_L2;0;Pgp_L1];
        
        localv_p=((delta_t*f')*(t\Pgp));
        gDp(ind1p:ind2p)=gDp(ind1p:ind2p)+localv_p;
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

% CC.D en el nodo S-E, fórmula de cuadratura para el término Pg de p en:
    % frontera Este
        Pgp_L7=int_dir_simpson(i,j,i,j+1,tt,3);
    % frontera Sur
        Pgp_L3=int_dir_simpson(i-1,j,i,j,tt,3);
    Pgp=(1/2)*[-Pgp_L7;Pgp_L3];
    
localv_p=((delta_t*f')*(t\Pgp));
gDp(ind1p)=gDp(ind1p)+localv_p;

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
    
    % CC.D en el nodo Oeste, fórmula de cuadratura para el término Pg de p en:
        % frontera Oeste 1
            Pgp_L4=int_dir_simpson(i,j,i,j-1,tt,3);
        % frontera Oeste 2
            Pgp_L11=int_dir_simpson(i,j,i,j+1,tt,3);
        Pgp=(1/2)*[Pgp_L4;0;Pgp_L11];
        
    localv_p=((delta_t*f')*(t\Pgp));
    gDp(ind1p)=gDp(ind1p)+localv_p(1);
    gDp(ind2p)=gDp(ind2p)+localv_p(2);
    
    % Central nodes: No hay CCD
    
    % East nodes
    i=N+1;
    ind2u=(i-1+(j-2)*N)*2;
    ind4u=(i-1+(j-1)*N)*2;
    
    ind1p=ind2u/2;
    ind2p=ind4u/2;
    
%     vdim=6;
%     vdim2=vdim/2;

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

       % CC.D en el nodo Este, fórmula de cuadratura para el término Pg de p en:
        % frontera Este 1
            Pgp_L7=int_dir_simpson(i,j,i,j-1,tt,3);
        % frontera Este 2
            Pgp_L14=int_dir_simpson(i,j,i,j+1,tt,3);
        Pgp=(1/2)*[-Pgp_L7;-Pgp_L14;0];
    
       localv_p=((delta_t*f')*(t\Pgp));
       gDp(ind1p)=gDp(ind1p)+localv_p(1);
       gDp(ind2p)=gDp(ind2p)+localv_p(2);
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

% CC.D en el nodo N-W, fórmula de cuadratura para el término Pg de p en:
    % frontera Oeste
        Pgp_L18=int_dir_simpson(i,j,i,j-1,tt,3);
    % frontera Norte
        Pgp_L22=int_dir_simpson(i,j,i+1,j,tt,3);
    Pgp=(1/2)*[Pgp_L18;-Pgp_L22];
    
localv_p=((delta_t*f')*(t\Pgp));
gDp(ind1p)=gDp(ind1p)+localv_p;

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
    
    % CC.D en el nodo Norte, fórmula de cuadratura para el término Pg de p en:
        % frontera Norte 1
            Pgp_L23=int_dir_simpson(i,j,i+1,j,tt,3);
        % frontera Norte 2
            Pgp_L22=int_dir_simpson(i-1,j,i,j,tt,3);
        Pgp=(1/2)*[0;-Pgp_L23;-Pgp_L22];
    
    localv_p=((delta_t*f')*(t\Pgp));
    gDp(ind1p:ind2p)=gDp(ind1p:ind2p)+localv_p;
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

% CC.D en el nodo N-E, fórmula de cuadratura para el término Pg de p en:
    % frontera Este
        Pgp_L21=int_dir_simpson(i,j,i,j-1,tt,3);
    % frontera Norte
        Pgp_L24=int_dir_simpson(i-1,j,i,j,tt,3);
    Pgp=(1/2)*[-Pgp_L21;-Pgp_L24];
    
localv_p=((delta_t*f')*(t\Pgp));
gDp(ind1p)=gDp(ind1p)+localv_p;
return
end
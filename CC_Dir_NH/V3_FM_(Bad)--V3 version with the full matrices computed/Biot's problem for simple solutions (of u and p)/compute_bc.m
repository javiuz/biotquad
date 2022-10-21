function [Pgu,Pgp]=compute_bc(tt)

global NN 

N=NN;

Pgu=zeros(8*N*(N+1),1);
Pgp=zeros(4*N*(N+1),1);

% South-West corner node
i=1;
j=1;

vdim=4;
vdim2=vdim/2;

ind_u=1:vdim;
ind_p=1:vdim2;

% CC.D en el nodo S-W, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
    % frontera Sur
        Pg1u_L1=int_dir_simpson(i,j,i+1,j,tt,1);
        Pg2u_L1=int_dir_simpson(i,j,i+1,j,tt,2);
    % frontera Oeste
        Pg1u_L4=int_dir_simpson(i,j,i,j+1,tt,1);
        Pg2u_L4=int_dir_simpson(i,j,i,j+1,tt,2);
        
    local_Pgu=[-Pg1u_L1;-Pg2u_L1;-Pg1u_L4;-Pg2u_L4];  
    Pgu(ind_u)=Pgu(ind_u) + local_Pgu;
   
% CC.D en el nodo S-W, fórmula de cuadratura para el término Pg de p en:
    % frontera Sur
        Pgp_L1=int_dir_simpson(i,j,i+1,j,tt,3);
    % frontera Oeste
        Pgp_L4=int_dir_simpson(i,j,i,j+1,tt,3);
    
        local_Pgp=[Pgp_L1;Pgp_L4];  
        Pgp(ind_p)=Pgp(ind_p) + local_Pgp;
        
ld_u=vdim;
ld_p=vdim2;

% South nodes (j=1)
vdim=6;
vdim2=vdim/2;

for i=2:N
     ind_u=ld_u+1:ld_u+vdim;
     ind_p=ld_p+1:ld_p+vdim2;
    
    % CC.D en el nodo Sur, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
        % frontera Sur 1
            Pg1u_L2=int_dir_simpson(i,j,i+1,j,tt,1);
            Pg2u_L2=int_dir_simpson(i,j,i+1,j,tt,2);
        % frontera Sur 2
            Pg1u_L1=int_dir_simpson(i-1,j,i,j,tt,1);
            Pg2u_L1=int_dir_simpson(i-1,j,i,j,tt,2);
        local_Pgu=[-Pg1u_L2;-Pg2u_L2;0;0;-Pg1u_L1;-Pg2u_L1];
    
        Pgu(ind_u)=Pgu(ind_u)+local_Pgu;
        
    % CC.D en los nodos Sur, fórmula de cuadratura para el término Pg de p en:
        % frontera Sur 1
            Pgp_L2=int_dir_simpson(i,j,i+1,j,tt,3); 
        % frontera Sur 2
            Pgp_L1=int_dir_simpson(i-1,j,i,j,tt,3); 
        local_Pgp=[Pgp_L2;0;Pgp_L1];
        
        Pgp(ind_p)=Pgp(ind_p)+local_Pgp;
        
        ld_u=ld_u+vdim;
        ld_p=ld_p+vdim2;
end

% South-East corner node (j=1)
i=N+1;

vdim=4;
vdim2=vdim/2;

ind_u=ld_u+1:ld_u+vdim;
ind_p=ld_p+1:ld_p+vdim2;

% CC.D en el nodo S-E, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
    % frontera Este
        Pg1u_L7=int_dir_simpson(i,j,i,j+1,tt,1);
        Pg2u_L7=int_dir_simpson(i,j,i,j+1,tt,2);
    % frontera Sur
        Pg1u_L3=int_dir_simpson(i-1,j,i,j,tt,1);
        Pg2u_L3=int_dir_simpson(i-1,j,i,j,tt,2);
    local_Pgu=[Pg1u_L7;Pg2u_L7;-Pg1u_L3;-Pg2u_L3];
    
Pgu(ind_u)=Pgu(ind_u)+local_Pgu;

% CC.D en el nodo S-E, fórmula de cuadratura para el término Pg de p en:
    % frontera Este
        Pgp_L7=int_dir_simpson(i,j,i,j+1,tt,3);
    % frontera Sur
        Pgp_L3=int_dir_simpson(i-1,j,i,j,tt,3);
    local_Pgp=[-Pgp_L7;Pgp_L3];
    
Pgp(ind_p)=Pgp(ind_p)+local_Pgp;

ld_u=ld_u+vdim;
ld_p=ld_p+vdim2;

for j=2:N

    % West nodes
    i=1;
    
    vdim=6;
    vdim2=vdim/2;
    
    ind_u=ld_u+1:ld_u+vdim;
    ind_p=ld_p+1:ld_p+vdim2;
    
    % CC.D en el nodo Oeste, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
        % frontera Oeste 1
            Pg1u_L4=int_dir_simpson(i,j,i,j-1,tt,1);
            Pg2u_L4=int_dir_simpson(i,j,i,j-1,tt,2);
        % frontera Oeste 2
            Pg1u_L11=int_dir_simpson(i,j,i,j+1,tt,1);
            Pg2u_L11=int_dir_simpson(i,j,i,j+1,tt,2);
        local_Pgu=[-Pg1u_L4;-Pg2u_L4;0;0;-Pg1u_L11;-Pg2u_L11];

    Pgu(ind_u)=Pgu(ind_u)+local_Pgu;
    
    % CC.D en el nodo Oeste, fórmula de cuadratura para el término Pg de p en:
        % frontera Oeste 1
            Pgp_L4=int_dir_simpson(i,j,i,j-1,tt,3);
        % frontera Oeste 2
            Pgp_L11=int_dir_simpson(i,j,i,j+1,tt,3);
        local_Pgp=[Pgp_L4;0;Pgp_L11];

    Pgp(ind_p)=Pgp(ind_p)+local_Pgp;
    
    ld_u=ld_u+vdim;
    ld_p=ld_p+vdim2;
    
    % Central nodes: 
% No hay CCD, pero tenemos que seguir avanzando en los índices para el 
% resto de C.C.
    vdim=8;
    vdim2=vdim/2;
    
    for i=2:N
        ld_u=ld_u+vdim;
        ld_p=ld_p+vdim2;
    end
    
    % East nodes
    i=N+1;
    
    vdim=6;
    vdim2=vdim/2;
    
    ind_u=ld_u+1:ld_u+vdim;
    ind_p=ld_p+1:ld_p+vdim2;

    % CC.D en el nodo Este, fórmula de cuadratura para los término Pg1 y Pg2 de u en:
        % frontera Este 1
            Pg1u_L7=int_dir_simpson(i,j,i,j-1,tt,1);
            Pg2u_L7=int_dir_simpson(i,j,i,j-1,tt,2);
        % frontera Este 2
            Pg1u_L14=int_dir_simpson(i,j,i,j+1,tt,1);
            Pg2u_L14=int_dir_simpson(i,j,i,j+1,tt,2);
        local_Pgu=[Pg1u_L7;Pg2u_L7;Pg1u_L14;Pg2u_L14;0;0];
        
       Pgu(ind_u)=Pgu(ind_u)+local_Pgu;
       
       % CC.D en el nodo Este, fórmula de cuadratura para el término Pg de p en:
        % frontera Este 1
            Pgp_L7=int_dir_simpson(i,j,i,j-1,tt,3);
        % frontera Este 2
            Pgp_L14=int_dir_simpson(i,j,i,j+1,tt,3);
        local_Pgp=[-Pgp_L7;-Pgp_L14;0];

       Pgp(ind_p)=Pgp(ind_p)+local_Pgp;
       
       ld_u=ld_u+vdim;
       ld_p=ld_p+vdim2;
end

% North-West corner node
i=1;
j=N+1;

vdim=4;
vdim2=vdim/2;

ind_u=ld_u+1:ld_u+vdim;
ind_p=ld_p+1:ld_p+vdim2;

% CC.D en el nodo N-W, fórmula de cuadratura para los término Pg1 y Pg2 en:
    % frontera Oeste
        Pg1u_L18=int_dir_simpson(i,j,i,j-1,tt,1);
        Pg2u_L18=int_dir_simpson(i,j,i,j-1,tt,2);
    % frontera Norte
        Pg1u_L22=int_dir_simpson(i,j,i+1,j,tt,1);
        Pg2u_L22=int_dir_simpson(i,j,i+1,j,tt,2);
    local_Pgu=[-Pg1u_L18;-Pg2u_L18;Pg1u_L22;Pg2u_L22];

Pgu(ind_u)=Pgu(ind_u)+local_Pgu;

% CC.D en el nodo N-W, fórmula de cuadratura para el término Pg de p en:
    % frontera Oeste
        Pgp_L18=int_dir_simpson(i,j,i,j-1,tt,3);
    % frontera Norte
        Pgp_L22=int_dir_simpson(i,j,i+1,j,tt,3);
    local_Pgp=[Pgp_L18;-Pgp_L22];
    
Pgp(ind_p)=Pgp(ind_p)+local_Pgp;

ld_u=ld_u+vdim;
ld_p=ld_p+vdim2;

% North nodes (j=N+1)
vdim=6;
vdim2=vdim/2;

for i=2:N
    ind_u=ld_u+1:ld_u+vdim;
    ind_p=ld_p+1:ld_p+vdim2;
    
    % CC.D en el nodo Norte, fórmula de cuadratura para los términos Pg1 y Pg2 de u en:
        % frontera Norte 1
            Pg1u_L23=int_dir_simpson(i,j,i+1,j,tt,1);
            Pg2u_L23=int_dir_simpson(i,j,i+1,j,tt,2);
        % frontera Norte 2
            Pg1u_L22=int_dir_simpson(i-1,j,i,j,tt,1);
            Pg2u_L22=int_dir_simpson(i-1,j,i,j,tt,2);
        local_Pgu=[0;0;Pg1u_L23;Pg2u_L23;Pg1u_L22;Pg2u_L22];
        
    Pgu(ind_u)=Pgu(ind_u)+local_Pgu;
    
    % CC.D en el nodo Norte, fórmula de cuadratura para el término Pg de p en:
        % frontera Norte 1
            Pgp_L23=int_dir_simpson(i,j,i+1,j,tt,3);
        % frontera Norte 2
            Pgp_L22=int_dir_simpson(i-1,j,i,j,tt,3);
        local_Pgp=[0;-Pgp_L23;-Pgp_L22];
        
    Pgp(ind_p)=Pgp(ind_p)+local_Pgp;
    
    ld_u=ld_u+vdim;
    ld_p=ld_p+vdim2;
end

% North-East corner node (j=N+1)
i=N+1;

vdim=4;
vdim2=vdim/2;

ind_u=ld_u+1:ld_u+vdim;
ind_p=ld_p+1:ld_p+vdim2;
          
% CC.D en el nodo N-E, fórmula de cuadratura para los términos Pg1 y Pg2 de u en:
    % frontera Este
        Pg1u_L21=int_dir_simpson(i,j,i,j-1,tt,1);
        Pg2u_L21=int_dir_simpson(i,j,i,j-1,tt,2);
    % frontera Norte
        Pg1u_L24=int_dir_simpson(i-1,j,i,j,tt,1);
        Pg2u_L24=int_dir_simpson(i-1,j,i,j,tt,2);
    local_Pgu=[Pg1u_L21;Pg2u_L21;Pg1u_L24;Pg2u_L24];

Pgu(ind_u)=Pgu(ind_u)+local_Pgu;

% CC.D en el nodo N-E, fórmula de cuadratura para el término Pg de p en:
    % frontera Este
        Pgp_L21=int_dir_simpson(i,j,i,j-1,tt,3);
    % frontera Norte
        Pgp_L24=int_dir_simpson(i-1,j,i,j,tt,3);
    local_Pgp=[-Pgp_L21;-Pgp_L24];

Pgp(ind_p)=Pgp(ind_p)+local_Pgp;

Pgu=(1/2)*Pgu;
Pgp=(1/2)*Pgp;
return
end
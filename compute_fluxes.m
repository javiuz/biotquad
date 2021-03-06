function [z,zx,zy]=compute_fluxes(paux,tt)

global NN 

N=NN;

z=zeros(4*N*(N+1),1);

% p=zeros(nx,ny);
% 
% for j=1:ny
%     for i=1:nx
%         p(i,j)=paux(i+(j-1)*ny);
%     end
% end

p=reshape(paux,N,N);

zx=zeros(N+1,N+1,2);
zy=zx;

% South-West corner node
i=1;
j=1;
vdim=2;
pdim=1;
t=zeros(vdim,vdim);
pp=zeros(pdim,1);
ld=1;

f=(1/2)*[1;1];

t(1,1)=kinv(2,2,i,j,i,j);
t(1,2)=kinv(2,1,i,j,i,j);
t(2,1)=kinv(1,2,i,j,i,j);
t(2,2)=kinv(1,1,i,j,i,j);
t=(1/4)*t;

% CC.D en el nodo S-W, fórmula de cuadratura para el término Pg de p en:
    % frontera Sur
        Pgp_L1=int_dir_simpson(i,j,i+1,j,tt,3);
    % frontera Oeste
        Pgp_L4=int_dir_simpson(i,j,i,j+1,tt,3);
    Pgp=(1/2)*[Pgp_L1;Pgp_L4];

pp(1,1)=p(i,j);
zz=t\(Pgp-f*pp);

z(ld:vdim)=zz;
ld=vdim;

zy(i,j,1)=zz(1,1);
zx(i,j,2)=zz(2,1);

% South nodes (j=1)
vdim=3;
pdim=2;
t=zeros(vdim,vdim);
pp=zeros(pdim,1);

% Matrix f is the same for every South node
f=(1/2)*[0 1;-1 1;1 0];

for i=2:N

    % if CC==D
    t(1,1)=kinv(2,2,i,j,i,j);
    t(1,2)=kinv(2,1,i,j,i,j);
%    t(1,3)=0;
    t(2,1)=kinv(1,2,i,j,i,j);
    t(2,2)=kinv(1,1,i-1,j,i,j)+kinv(1,1,i,j,i,j);
    t(2,3)=kinv(1,2,i-1,j,i,j);
%    t(3,1)=0;
    t(3,2)=kinv(2,1,i-1,j,i,j);
    t(3,3)=kinv(2,2,i-1,j,i,j);
    t=(1/4)*t;
    
    % CC.D en los nodos Sur, fórmula de cuadratura para el término Pg de p en:
        % frontera Sur 1
            Pgp_L2=int_dir_simpson(i,j,i+1,j,tt,3); 
        % frontera Sur 2
            Pgp_L1=int_dir_simpson(i-1,j,i,j,tt,3); 
        Pgp=(1/2)*[Pgp_L2;0;Pgp_L1];
    
    pp(1,1)=p(i-1,j);
    pp(2,1)=p(i,j);
    zz=t\(Pgp-f*pp);
    
    z(ld+1:ld+vdim)=zz;
    ld=ld+vdim;

    zy(i,j,1)=zz(1,1);
    zx(i,j,2)=zz(2,1);
    zy(i,j,2)=zz(3,1);
end 
    
% South-East corner node (j=1)
i=N+1;
vdim=2;
pdim=1;
t=zeros(vdim,vdim);
pp=zeros(pdim,1);

f=(1/2)*[-1;1];

t(1,1)=kinv(1,1,i-1,j,i,j);
t(1,2)=kinv(1,2,i-1,j,i,j);
t(2,1)=kinv(2,1,i-1,j,i,j);
t(2,2)=kinv(2,2,i-1,j,i,j);
t=(1/4)*t;

% CC.D en el nodo S-E, fórmula de cuadratura para el término Pg de p en:
    % frontera Este
        Pgp_L7=int_dir_simpson(i,j,i,j+1,tt,3);
    % frontera Sur
        Pgp_L3=int_dir_simpson(i-1,j,i,j,tt,3);
    Pgp=(1/2)*[-Pgp_L7;Pgp_L3];

pp(1,1)=p(i-1,j);
zz=t\(Pgp-f*pp);

z(ld+1:ld+vdim)=zz;
ld=ld+vdim;

zx(i,j,2)=zz(1,1);
zy(i,j,2)=zz(2,1);

for j=2:N

    % West nodes
    i=1;
    vdim=3;
    pdim=2;
    t=zeros(vdim,vdim);
    pp=zeros(pdim,1);
    
    f=(1/2)*[1 0;-1 1;0 1];
    
    t(1,1)=kinv(1,1,i,j-1,i,j);
    t(1,2)=kinv(1,2,i,j-1,i,j);
%    t(1,3)=0.
    t(2,1)=kinv(2,1,i,j-1,i,j);
    t(2,2)=kinv(2,2,i,j-1,i,j)+kinv(2,2,i,j,i,j);
    t(2,3)=kinv(2,1,i,j,i,j);
%    t(3,1)=0.
    t(3,2)=kinv(1,2,i,j,i,j);
    t(3,3)=kinv(1,1,i,j,i,j);
    t=(1/4)*t;
    
    % CC.D en el nodo Oeste, fórmula de cuadratura para el término Pg de p en:
        % frontera Oeste 1
            Pgp_L4=int_dir_simpson(i,j,i,j-1,tt,3);
        % frontera Oeste 2
            Pgp_L11=int_dir_simpson(i,j,i,j+1,tt,3);
        Pgp=(1/2)*[Pgp_L4;0;Pgp_L11];
        
    pp(1,1)=p(i,j-1);
    pp(2,1)=p(i,j);
    zz=t\(Pgp-f*pp);
    
    z(ld+1:ld+vdim)=zz;
    ld=ld+vdim;
    
    zx(i,j,1)=zz(1,1);
    zy(i,j,1)=zz(2,1);
    zx(i,j,2)=zz(3,1);

    % Central nodes 
    vdim=4;
    pdim=4;
    t=zeros(vdim,vdim);
    pp=zeros(pdim,1);

    % Matrix f is the same for every central node
    f=(1/2)*[-1 1 0 0;0 -1 0 1;0 0 -1 1;-1 0 1 0];
    
    for i=2:N

        t(1,1)=kinv(1,1,i-1,j-1,i,j)+kinv(1,1,i,j-1,i,j);
        t(1,2)=kinv(1,2,i,j-1,i,j);
        %t(1,3)=0._8
        t(1,4)=kinv(1,2,i-1,j-1,i,j);
        t(2,1)=kinv(2,1,i,j-1,i,j);
        t(2,2)=kinv(2,2,i,j-1,i,j)+kinv(2,2,i,j,i,j);
        t(2,3)=kinv(2,1,i,j,i,j);
        %t(2,4)=0._8
        %t(3,1)=0._8
        t(3,2)=kinv(1,2,i,j,i,j);
        t(3,3)=kinv(1,1,i,j,i,j)+kinv(1,1,i-1,j,i,j);
        t(3,4)=kinv(1,2,i-1,j,i,j);
        t(4,1)=kinv(2,1,i-1,j-1,i,j);
        %t(4,2)=0._8
        t(4,3)=kinv(2,1,i-1,j,i,j);
        t(4,4)=kinv(2,2,i-1,j,i,j)+kinv(2,2,i-1,j-1,i,j);
        t=(1/4)*t;
        
        pp(1,1)=p(i-1,j-1);
        pp(2,1)=p(i,j-1);
        pp(3,1)=p(i-1,j);
        pp(4,1)=p(i,j);
        zz=-t\(f*pp);
        
        z(ld+1:ld+vdim)=zz;
        ld=ld+vdim;

        zx(i,j,1)=zz(1,1);
        zy(i,j,1)=zz(2,1);
        zx(i,j,2)=zz(3,1);
        zy(i,j,2)=zz(4,1);
    end 

    % East nodes
    i=N+1;

    vdim=3;
    pdim=2;
    t=zeros(vdim,vdim);
    pp=zeros(pdim,1);
    f=(1/2)*[-1 0;0 -1;-1 1];

    t(1,1)=kinv(1,1,i-1,j-1,i,j);
    %t(1,2)=0._8
    t(1,3)=kinv(1,2,i-1,j-1,i,j);
    %t(2,1)=0._8
    t(2,2)=kinv(1,1,i-1,j,i,j);
    t(2,3)=kinv(1,2,i-1,j,i,j);
    t(3,1)=kinv(2,1,i-1,j-1,i,j);
    t(3,2)=kinv(2,1,i-1,j,i,j);
    t(3,3)=kinv(2,2,i-1,j-1,i,j)+kinv(2,2,i-1,j,i,j);
    t=(1/4)*t;
    
    % CC.D en el nodo Este, fórmula de cuadratura para el término Pg de p en:
        % frontera Este 1
            Pgp_L7=int_dir_simpson(i,j,i,j-1,tt,3);
        % frontera Este 2
            Pgp_L14=int_dir_simpson(i,j,i,j+1,tt,3);
        Pgp=(1/2)*[-Pgp_L7;-Pgp_L14;0];

    pp(1,1)=p(i-1,j-1);
    pp(2,1)=p(i-1,j);
    zz=t\(Pgp-f*pp);
    
    z(ld+1:ld+vdim)=zz;
    ld=ld+vdim;

    zx(i,j,1)=zz(1,1);
    zx(i,j,2)=zz(2,1);
    zy(i,j,2)=zz(3,1);
    
end

% North-West corner node
i=1;
j=N+1;
vdim=2;
pdim=1;
t=zeros(vdim,vdim);
pp=zeros(pdim,1);

f=(1/2)*[1;-1];

t(1,1)=kinv(1,1,i,j-1,i,j);
t(1,2)=kinv(1,2,i,j-1,i,j);
t(2,1)=kinv(2,1,i,j-1,i,j);
t(2,2)=kinv(2,2,i,j-1,i,j);
t=(1/4)*t;

% CC.D en el nodo N-W, fórmula de cuadratura para el término Pg de p en:
    % frontera Oeste
        Pgp_L18=int_dir_simpson(i,j,i,j-1,tt,3);
    % frontera Norte
        Pgp_L22=int_dir_simpson(i,j,i+1,j,tt,3);
    Pgp=(1/2)*[Pgp_L18;-Pgp_L22];

pp(1,1)=p(i,j-1);
zz=t\(Pgp-f*pp);

z(ld+1:ld+vdim)=zz;
ld=ld+vdim;

zx(i,j,1)=zz(1,1);
zy(i,j,1)=zz(2,1);

% North nodes (j=ny+1)
vdim=3;
pdim=2;
t=zeros(vdim,vdim);
pp=zeros(pdim,1);

% Matrix f is the same for every North node
f=(1/2)*[-1 1;0 -1;-1 0];

for i=2:N

    t(1,1)=kinv(1,1,i-1,j-1,i,j)+kinv(1,1,i,j-1,i,j);
    t(1,2)=kinv(1,2,i,j-1,i,j);
    t(1,3)=kinv(1,2,i-1,j-1,i,j);
    t(2,1)=kinv(2,1,i,j-1,i,j);
    t(2,2)=kinv(2,2,i,j-1,i,j);
    %t(2,3)=0._8
    t(3,1)=kinv(2,1,i-1,j-1,i,j);
    %t(3,2)=0._8
    t(3,3)=kinv(2,2,i-1,j-1,i,j);
    t=(1/4)*t;
    
    % CC.D en el nodo Norte, fórmula de cuadratura para el término Pg de p en:
        % frontera Norte 1
            Pgp_L23=int_dir_simpson(i,j,i+1,j,tt,3);
        % frontera Norte 2
            Pgp_L22=int_dir_simpson(i-1,j,i,j,tt,3);
        Pgp=(1/2)*[0;-Pgp_L23;-Pgp_L22];
  
    pp(1,1)=p(i-1,j-1);
    pp(2,1)=p(i,j-1);
    zz=t\(Pgp-f*pp);
    
    z(ld+1:ld+vdim)=zz;
    ld=ld+vdim;

    zx(i,j,1)=zz(1,1);
    zy(i,j,1)=zz(2,1);
    zy(i,j,2)=zz(3,1);
end 

% North-East corner node (j=ny+1)
i=N+1;
vdim=2;
pdim=1;
t=zeros(vdim,vdim);
pp=zeros(pdim,1);

f=(1/2)*[-1;-1];

t(1,1)=kinv(1,1,i-1,j-1,i,j);
t(1,2)=kinv(1,2,i-1,j-1,i,j);
t(2,1)=kinv(2,1,i-1,j-1,i,j);
t(2,2)=kinv(2,2,i-1,j-1,i,j);
t=(1/4)*t;

% CC.D en el nodo N-E, fórmula de cuadratura para el término Pg de p en:
    % frontera Este
        Pgp_L21=int_dir_simpson(i,j,i,j-1,tt,3);
    % frontera Norte
        Pgp_L24=int_dir_simpson(i-1,j,i,j,tt,3);
    Pgp=(1/2)*[-Pgp_L21;-Pgp_L24];

pp(1,1)=p(i-1,j-1);
zz=t\(Pgp-f*pp);

z(ld+1:ld+vdim)=zz;

zx(i,j,1)=zz(1,1);
zy(i,j,2)=zz(2,1);

return
end
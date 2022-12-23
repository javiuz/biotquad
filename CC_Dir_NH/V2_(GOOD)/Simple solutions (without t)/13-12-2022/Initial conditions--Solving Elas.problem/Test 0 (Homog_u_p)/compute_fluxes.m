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

pp(1,1)=p(i,j);
zz=t\(-f*pp);

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
    
    pp(1,1)=p(i-1,j);
    pp(2,1)=p(i,j);
    zz=t\(-f*pp);
    
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

pp(1,1)=p(i-1,j);
zz=t\(-f*pp);

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
        
    pp(1,1)=p(i,j-1);
    pp(2,1)=p(i,j);
    zz=t\(-f*pp);
    
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

    pp(1,1)=p(i-1,j-1);
    pp(2,1)=p(i-1,j);
    zz=t\(-f*pp);
    
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

pp(1,1)=p(i,j-1);
zz=t\(-f*pp);

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
  
    pp(1,1)=p(i-1,j-1);
    pp(2,1)=p(i,j-1);
    zz=t\(-f*pp);
    
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

pp(1,1)=p(i-1,j-1);
zz=t\(-f*pp);

z(ld+1:ld+vdim)=zz;

zx(i,j,1)=zz(1,1);
zy(i,j,2)=zz(2,1);

return
end
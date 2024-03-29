function sigma=build_sigma_0_cartesian(tt,nu)

global NN x y 

N=NN;

sigma=zeros(8*N*(N+1),1);

% South-West corner node
i=1;
j=1;

vdim=4;

x1=x(i,j);  
y1=y(i,j);            
% x2=x(i+1,j);
% y2=y(i+1,j);
% x3=x(i+1,j+1);
% y3=y(i+1,j+1);
% x4=x(i,j+1);
% y4=y(i,j+1);

ss=zeros(vdim,1);
ss(1)=sol_exactax_sigma(x1,y1,tt,nu,1,2);
ss(2)=sol_exactax_sigma(x1,y1,tt,nu,2,2);
ss(3)=sol_exactax_sigma(x1,y1,tt,nu,1,1);
ss(4)=sol_exactax_sigma(x1,y1,tt,nu,2,1);

indr=1:vdim;
sigma(indr)=ss;
ld=vdim;

% South nodes (j=1)
vdim=6;
ss=zeros(vdim,1);

for i=2:N
   
%     x1=x(i-1,j);
%     y1=y(i-1,j);
    x2=x(i,j);
    y2=y(i,j);
%     x3=x(i,j+1);
%     y3=y(i,j+1);
%     x4=x(i-1,j+1);
%     y4=y(i-1,j+1);
%     x5=x(i+1,j);
%     y5=y(i+1,j);
%     x6=x(i+1,j+1);
%     y6=y(i+1,j+1);
    
    ss(1)=sol_exactax_sigma(x2,y2,tt,nu,1,2);
    ss(2)=sol_exactax_sigma(x2,y2,tt,nu,2,2);
    ss(3)=sol_exactax_sigma(x2,y2,tt,nu,1,1);
    ss(4)=sol_exactax_sigma(x2,y2,tt,nu,2,1);
    ss(5)=ss(1);
    ss(6)=ss(2);
    
    indr=ld+1:ld+vdim;
    sigma(indr)=ss;
    ld=ld+vdim;
end

% South-East corner node (j=1)
i=N+1;
vdim=4;
ss=zeros(vdim,1);

% x1=x(i-1,j);  
% y1=y(i-1,j);            
x2=x(i,j);
y2=y(i,j);
% x3=x(i,j+1);
% y3=y(i,j+1);
% x4=x(i-1,j+1);
% y4=y(i-1,j+1);

ss(1)=sol_exactax_sigma(x2,y2,tt,nu,1,1);
ss(2)=sol_exactax_sigma(x2,y2,tt,nu,2,1);
ss(3)=sol_exactax_sigma(x2,y2,tt,nu,1,2);
ss(4)=sol_exactax_sigma(x2,y2,tt,nu,2,2);

indr=ld+1:ld+vdim;
sigma(indr)=ss;
ld=ld+vdim;

for j=2:N

    % West nodes
    i=1;
    
    vdim=6;
    ss=zeros(vdim,1);
    
%     x1=x(i,j-1);
%     y1=y(i,j-1);
%     x2=x(i+1,j-1);
%     y2=y(i+1,j-1);
%     x3=x(i+1,j);
%     y3=y(i+1,j);
    x4=x(i,j);
    y4=y(i,j);
%     x5=x(i+1,j+1);
%     y5=y(i+1,j+1);
%     x6=x(i,j+1);
%     y6=y(i,j+1);
    
    ss(1)=sol_exactax_sigma(x4,y4,tt,nu,1,1);
    ss(2)=sol_exactax_sigma(x4,y4,tt,nu,2,1);
    ss(3)=sol_exactax_sigma(x4,y4,tt,nu,1,2);
    ss(4)=sol_exactax_sigma(x4,y4,tt,nu,2,2);
    ss(5)=ss(1);
    ss(6)=ss(2);

    indr=ld+1:ld+vdim;
    sigma(indr)=ss;
    ld=ld+vdim;
    
    % Central nodes 
    vdim=8;

    ss=zeros(vdim,1);
         
    for i=2:N
        
%     x1=x(i-1,j-1);
%     y1=y(i-1,j-1);
%     x2=x(i,j-1);
%     y2=y(i,j-1);
    x3=x(i,j);
    y3=y(i,j);
%     x4=x(i-1,j);
%     y4=y(i-1,j);
%     x5=x(i+1,j-1);
%     y5=y(i+1,j-1);
%     x6=x(i+1,j);
%     y6=y(i+1,j);
%     x7=x(i+1,j+1);
%     y7=y(i+1,j+1);
%     x8=x(i,j+1);
%     y8=y(i,j+1);
%     x9=x(i-1,j+1);
%     y9=y(i-1,j+1);
    
    ss(1)=sol_exactax_sigma(x3,y3,tt,nu,1,1);
    ss(2)=sol_exactax_sigma(x3,y3,tt,nu,2,1);
    ss(3)=sol_exactax_sigma(x3,y3,tt,nu,1,2);
    ss(4)=sol_exactax_sigma(x3,y3,tt,nu,2,2);
    ss(5)=ss(1);
    ss(6)=ss(2);
    ss(7)=ss(3);
    ss(8)=ss(4);

    indr=ld+1:ld+vdim;
    sigma(indr)=ss;
    ld=ld+vdim;
    end
    
    % East nodes
    i=N+1;
    
    vdim=6;
    ss=zeros(vdim,1);
    
%     x1=x(i-1,j-1);
%     y1=y(i-1,j-1);
%     x2=x(i,j-1);
%     y2=y(i,j-1);
    x3=x(i,j);
    y3=y(i,j);
%     x4=x(i-1,j);
%     y4=y(i-1,j);
%     x5=x(i,j+1);
%     y5=y(i,j+1);
%     x6=x(i-1,j+1);
%     y6=y(i-1,j+1);   
    
    ss(1)=sol_exactax_sigma(x3,y3,tt,nu,1,1);
    ss(2)=sol_exactax_sigma(x3,y3,tt,nu,2,1);
    ss(3)=ss(1);
    ss(4)=ss(2);
    ss(5)=sol_exactax_sigma(x3,y3,tt,nu,1,2);
    ss(6)=sol_exactax_sigma(x3,y3,tt,nu,2,2);

    indr=ld+1:ld+vdim;
    sigma(indr)=ss;
    ld=ld+vdim;
end

% North-West corner node
i=1;
vdim=4;
ss=zeros(vdim,1);

% x1=x(i,j-1);
% y1=y(i,j-1);           
% x2=x(i+1,j-1);
% y2=y(i+1,j-1);
% x3=x(i+1,j);
% y3=y(i+1,j);
x4=x(i,j);
y4=y(i,j);

ss(1)=sol_exactax_sigma(x4,y4,tt,nu,1,1);
ss(2)=sol_exactax_sigma(x4,y4,tt,nu,2,1);
ss(3)=sol_exactax_sigma(x4,y4,tt,nu,1,2);
ss(4)=sol_exactax_sigma(x4,y4,tt,nu,2,2);

indr=ld+1:ld+vdim;
sigma(indr)=ss;
ld=ld+vdim;

% North nodes (j=N+1)
vdim=6;
ss=zeros(vdim,1);

for i=2:N
    
%     x1=x(i-1,j-1);
%     y1=y(i-1,j-1);
%     x2=x(i,j-1);
%     y2=y(i,j-1);
    x3=x(i,j);
    y3=y(i,j);
%     x4=x(i-1,j);
%     y4=y(i-1,j);
%     x5=x(i+1,j-1);
%     y5=y(i+1,j-1);
%     x6=x(i+1,j);
%     y6=y(i+1,j);   
    
    ss(1)=sol_exactax_sigma(x3,y3,tt,nu,1,1);
    ss(2)=sol_exactax_sigma(x3,y3,tt,nu,2,1);
    ss(3)=sol_exactax_sigma(x3,y3,tt,nu,1,2);
    ss(4)=sol_exactax_sigma(x3,y3,tt,nu,2,2);
    ss(5)=ss(3);
    ss(6)=ss(4);

    indr=ld+1:ld+vdim;
    sigma(indr)=ss;
    ld=ld+vdim;
end

% North-East corner node (j=N+1)
i=N+1;
vdim=4;
ss=zeros(vdim,1);

% x1=x(i-1,j-1);  
% y1=y(i-1,j-1);            
% x2=x(i,j-1);
% y2=y(i,j-1);
x3=x(i,j);
y3=y(i,j);
% x4=x(i-1,j);
% y4=y(i-1,j);

ss(1)=sol_exactax_sigma(x3,y3,tt,nu,1,1);
ss(2)=sol_exactax_sigma(x3,y3,tt,nu,2,1);
ss(3)=sol_exactax_sigma(x3,y3,tt,nu,1,2);
ss(4)=sol_exactax_sigma(x3,y3,tt,nu,2,2);

indr=ld+1:ld+vdim;
sigma(indr)=ss;
sigma=sigma*(1/N);
return
end
function sigma=build_exact_sol_sigma(tt)

global NN x y 

N=NN;

sigma=zeros(8*N*(N+1),1);

% South-West corner node
i=1;
j=1;
vdim=4;

xv=x(i,j);  
yv=y(i,j);         
    
sigma(1,1)=sol_exactax_sigma(xv,yv,tt,tt,1,2);
sigma(2,1)=sol_exactax_sigma(xv,yv,tt,tt,2,2);
sigma(3,1)=sol_exactax_sigma(xv,yv,tt,tt,1,1);
sigma(4,1)=sol_exactax_sigma(xv,yv,tt,tt,2,1);
ld=vdim;

% South nodes (j=1)
vdim=6;

for i=2:N
    xv=x(i,j);
    yv=y(i,j);
    
    sigma(ld+1,1)=sol_exactax_sigma(xv,yv,tt,tt,1,2);
    sigma(ld+2,1)=sol_exactax_sigma(xv,yv,tt,tt,2,2);
    sigma(ld+3,1)=sol_exactax_sigma(xv,yv,tt,tt,1,1);
    sigma(ld+4,1)=sol_exactax_sigma(xv,yv,tt,tt,2,1);
    sigma(ld+5,1)=sigma(ld+1,1);
    sigma(ld+6,1)=sigma(ld+2,1);
    ld=ld+vdim;
end

% South-East corner node (j=1)
i=N+1;
vdim=4;
       
xv=x(i,j);
yv=y(i,j);

sigma(ld+1,1)=sol_exactax_sigma(xv,yv,tt,tt,1,1);
sigma(ld+2,1)=sol_exactax_sigma(xv,yv,tt,tt,2,1);
sigma(ld+3,1)=sol_exactax_sigma(xv,yv,tt,tt,1,2);
sigma(ld+4,1)=sol_exactax_sigma(xv,yv,tt,tt,2,2);

ld=ld+vdim;

for j=2:N

    % West nodes
    i=1;
    
    vdim=6;
    
    xv=x(i,j);
    yv=y(i,j);
    
    sigma(ld+1,1)=sol_exactax_sigma(xv,yv,tt,tt,1,1);
    sigma(ld+2,1)=sol_exactax_sigma(xv,yv,tt,tt,2,1);
    sigma(ld+3,1)=sol_exactax_sigma(xv,yv,tt,tt,1,2);
    sigma(ld+4,1)=sol_exactax_sigma(xv,yv,tt,tt,2,2);
    sigma(ld+5,1)=sigma(ld+1,1);
    sigma(ld+6,1)=sigma(ld+2,1);
    
    ld=ld+vdim;
    
    % Central nodes 
    vdim=8;
         
    for i=2:N     
        xv=x(i,j);
        yv=y(i,j);
    
        sigma(ld+1)=sol_exactax_sigma(xv,yv,tt,tt,1,1);
        sigma(ld+2)=sol_exactax_sigma(xv,yv,tt,tt,2,1);
        sigma(ld+3)=sol_exactax_sigma(xv,yv,tt,tt,1,2);
        sigma(ld+4)=sol_exactax_sigma(xv,yv,tt,tt,2,2);
        sigma(ld+5)=sigma(ld+1);
        sigma(ld+6)=sigma(ld+2);
        sigma(ld+7)=sigma(ld+3);
        sigma(ld+8)=sigma(ld+4);
        
        ld=ld+vdim;
    end
    
    % East nodes
    i=N+1;
    
    vdim=6;
    
    xv=x(i,j);
    yv=y(i,j);
    
    sigma(ld+1,1)=sol_exactax_sigma(xv,yv,tt,tt,1,1);
    sigma(ld+2,1)=sol_exactax_sigma(xv,yv,tt,tt,2,1);
    sigma(ld+3,1)=sigma(ld+1,1);
    sigma(ld+4,1)=sigma(ld+2,1);
    sigma(ld+5,1)=sol_exactax_sigma(xv,yv,tt,tt,1,2);
    sigma(ld+6,1)=sol_exactax_sigma(xv,yv,tt,tt,2,2);
    
    ld=ld+vdim;
end

% North-West corner node
i=1;
j=N+1;

vdim=4;

xv=x(i,j);
yv=y(i,j);

sigma(ld+1,1)=sol_exactax_sigma(xv,yv,tt,tt,1,1);
sigma(ld+2,1)=sol_exactax_sigma(xv,yv,tt,tt,2,1);
sigma(ld+3,1)=sol_exactax_sigma(xv,yv,tt,tt,1,2);
sigma(ld+4,1)=sol_exactax_sigma(xv,yv,tt,tt,2,2);

ld=ld+vdim;

% North nodes (j=N+1)
vdim=6;

for i=2:N
    xv=x(i,j);
    yv=y(i,j);
    
    sigma(ld+1,1)=sol_exactax_sigma(xv,yv,tt,tt,1,1);
    sigma(ld+2,1)=sol_exactax_sigma(xv,yv,tt,tt,2,1);
    sigma(ld+3,1)=sol_exactax_sigma(xv,yv,tt,tt,1,2);
    sigma(ld+4,1)=sol_exactax_sigma(xv,yv,tt,tt,2,2);
    sigma(ld+5,1)=sigma(ld+3,1);
    sigma(ld+6,1)=sigma(ld+4,1);
    
    ld=ld+vdim;
end

% North-East corner node (j=N+1)
i=N+1;
% vdim=4;

xv=x(i,j);
yv=y(i,j);

sigma(ld+1,1)=sol_exactax_sigma(xv,yv,tt,tt,1,1);
sigma(ld+2,1)=sol_exactax_sigma(xv,yv,tt,tt,2,1);
sigma(ld+3,1)=sol_exactax_sigma(xv,yv,tt,tt,1,2);
sigma(ld+4,1)=sol_exactax_sigma(xv,yv,tt,tt,2,2);

return
end
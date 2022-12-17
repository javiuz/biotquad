function sigma=build_sigma_0_gral_grid(tt,nu)

global NN x y 

N=NN;

sigma=zeros(8*N*(N+1),1);

sx1=zeros(N+1,N+1,2);
sx2=sx1;
sy1=sx1;
sy2=sx1;

normal_n=sx1;
normal_s=normal_n;
normal_e=normal_n;
normal_w=normal_n;

% South-West corner node

% Comp. normales de las tensiones (1ª fila por un lado y 2ª fila por otro) 
% vienen dadas por normal_e y normal_n (celda noreste de los apuntes). 

i=1;
j=1;

% Celda NE

normal_n(i,j,1)=(y(i,j+1)-y(i,j))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
normal_n(i,j,2)=(x(i,j)-x(i,j+1))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
normal_e(i,j,1)=(y(i,j)-y(i+1,j))/sqrt((y(i,j)-y(i+1,j))^2+(x(i+1,j)-x(i,j))^2);
normal_e(i,j,2)=(x(i+1,j)-x(i,j))/sqrt((y(i,j)-y(i+1,j))^2+(x(i+1,j)-x(i,j))^2);

denom_ne=normal_e(i,j,1)*normal_n(i,j,2)-normal_n(i,j,1)*normal_e(i,j,2);
    alpha_ne=normal_n(i,j,2)/denom_ne;
    beta_ne=-normal_e(i,j,2)/denom_ne;
    gamma_ne=-normal_n(i,j,1)/denom_ne;
    delta_ne=normal_e(i,j,1)/denom_ne;

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

sy1(i,j,1)=sol_exactax_sigma(x1,y1,tt,nu,1,2);
sy2(i,j,1)=sol_exactax_sigma(x1,y1,tt,nu,2,2);
sx1(i,j,2)=sol_exactax_sigma(x1,y1,tt,nu,1,1);
sx2(i,j,2)=sol_exactax_sigma(x1,y1,tt,nu,2,1);

    % 1ª fila de la matriz de tensiones
        % Comp. horizontal del flujo 
    s_flux_h_ne_1=alpha*sy1(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                beta*sx1(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
        % Comp. vertical del flujo 
    s_flux_v_ne_1=gamma*sy1(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                delta*sx1(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
    % 2ª fila de la matriz de tensiones
        % Comp. horizontal del flujo         
    s_flux_h_ne_2=alpha*sy2(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                beta*sx2(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
        % Comp. vertical del flujo 
    s_flux_v_ne_2=gamma*sy2(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                delta*sx2(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
        
ss(1)=s_flux_v_ne_1;
ss(2)=s_flux_v_ne_2;
ss(3)=s_flux_h_ne_1;
ss(4)=s_flux_h_ne_2;
            
% ss(1)=sy1(i,j,1);
% ss(2)=sy2(i,j,1);
% ss(3)=sx1(i,j,2);
% ss(4)=sx2(i,j,2);

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

    normal_n(i,j,1)=(y(i,j+1)-y(i,j))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
    normal_n(i,j,2)=(x(i,j)-x(i,j+1))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
    normal_e(i,j,1)=(y(i,j)-y(i+1,j))/sqrt((y(i,j)-y(i+1,j))^2+(x(i+1,j)-x(i,j))^2);
    normal_e(i,j,2)=(x(i+1,j)-x(i,j))/sqrt((y(i,j)-y(i+1,j))^2+(x(i+1,j)-x(i,j))^2);
    normal_w(i,j,1)=(y(i-1,j)-y(i,j))/sqrt((y(i-1,j)-y(i,j))^2+(x(i,j)-x(i-1,j))^2);
    normal_w(i,j,2)=(x(i,j)-x(i-1,j))/sqrt((y(i-1,j)-y(i,j))^2+(x(i,j)-x(i-1,j))^2);
    
    denom_ne=normal_e(i,j,1)*normal_n(i,j,2)-normal_n(i,j,1)*normal_e(i,j,2);
        alpha_ne=normal_n(i,j,2)/denom_ne;
        beta_ne=-normal_e(i,j,2)/denom_ne;
        gamma_ne=-normal_n(i,j,1)/denom_ne;
        delta_ne=normal_e(i,j,1)/denom_ne;
        
    denom_nw=normal_n(i,j,1)*normal_w(i,j,2)-normal_w(i,j,1)*normal_n(i,j,2);
        alpha_nw=normal_w(i,j,2)/denom_nw;
        beta_nw=-normal_n(i,j,2)/denom_nw;
        gamma_nw=-normal_w(i,j,1)/denom_nw;
        delta_nw=normal_n(i,j,1)/denom_nw;
    
    sy1(i,j,1)=sol_exactax_sigma(x2,y2,tt,nu,1,2);
    sy2(i,j,1)=sol_exactax_sigma(x2,y2,tt,nu,2,2);
    sx1(i,j,2)=sol_exactax_sigma(x2,y2,tt,nu,1,1);
    sx2(i,j,2)=sol_exactax_sigma(x2,y2,tt,nu,2,1);
    sy1(i,j,2)=sy1(i,j,1);
    sy2(i,j,2)=sy2(i,j,1);
    
    % Celda NE
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_ne_1=alpha_ne*sy1(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                   beta_ne*sx1(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
            % Comp. vertical del flujo 
     s_flux_v_ne_1=gamma_ne*sy1(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                   delta_ne*sx1(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);   
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_ne_2=alpha_ne*sy2(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                   beta_ne*sx2(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
            % Comp. vertical del flujo 
     s_flux_v_ne_2=gamma_ne*sy2(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                   delta_ne*sx2(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);   

    % Celda NO
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_nw_1=alpha_nw*sx1(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   beta_nw*sy1(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_nw_1=gamma_nw*sx1(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   delta_nw*sy1(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);   
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_nw_2=alpha_nw*sx2(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   beta_nw*sy2(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_nw_2=gamma_nw*sx2(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   delta_nw*sy2(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);   
               
    ss(1)=s_flux_v_ne_1;
    ss(2)=s_flux_v_ne_2;
    ss(3)=(s_flux_h_ne_1 + s_flux_h_nw_1)/2;
    ss(4)=(s_flux_h_ne_2 + s_flux_h_nw_2)/2;
    ss(5)=s_flux_v_nw_1;
    ss(6)=s_flux_v_nw_2;
               
%     ss(1)=sol_exactax_sigma(x2,y2,tt,nu,1,2);
%     ss(2)=sol_exactax_sigma(x2,y2,tt,nu,2,2);
%     ss(3)=sol_exactax_sigma(x2,y2,tt,nu,1,1);
%     ss(4)=sol_exactax_sigma(x2,y2,tt,nu,2,1);
%     ss(5)=ss(1);
%     ss(6)=ss(2);
    
    indr=ld+1:ld+vdim;
    sigma(indr)=ss;
    ld=ld+vdim;
end

% South-East corner node (j=1)
i=N+1;
vdim=4;
ss=zeros(vdim,1);

% Celda NO
normal_n(i,j,1)=(y(i,j+1)-y(i,j))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
normal_n(i,j,2)=(x(i,j)-x(i,j+1))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
normal_w(i,j,1)=(y(i-1,j)-y(i,j))/sqrt((y(i-1,j)-y(i,j))^2+(x(i,j)-x(i-1,j))^2);
normal_w(i,j,2)=(x(i,j)-x(i-1,j))/sqrt((y(i-1,j)-y(i,j))^2+(x(i,j)-x(i-1,j))^2);

denom_nw=normal_n(i,j,1)*normal_w(i,j,2)-normal_w(i,j,1)*normal_n(i,j,2);
    alpha_nw=normal_w(i,j,2)/denom_nw;
    beta_nw=-normal_n(i,j,2)/denom_nw;
    gamma_nw=-normal_w(i,j,1)/denom_nw;
    delta_nw=normal_n(i,j,1)/denom_nw;

% x1=x(i-1,j);  
% y1=y(i-1,j);            
x2=x(i,j);
y2=y(i,j);
% x3=x(i,j+1);
% y3=y(i,j+1);
% x4=x(i-1,j+1);
% y4=y(i-1,j+1);

sx1(i,j,2)=sol_exactax_sigma(x2,y2,tt,nu,1,1);
sx2(i,j,2)=sol_exactax_sigma(x2,y2,tt,nu,2,1);
sy1(i,j,2)=sol_exactax_sigma(x2,y2,tt,nu,1,2);
sy2(i,j,2)=sol_exactax_sigma(x2,y2,tt,nu,2,2);

    % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_nw_1=alpha_nw*sx1(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   beta_nw*sy1(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_nw_1=gamma_nw*sx1(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   delta_nw*sy1(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);   
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_nw_2=alpha_nw*sx2(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   beta_nw*sy2(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_nw_2=gamma_nw*sx2(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   delta_nw*sy2(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);  

ss(1)=s_flux_h_nw_1;
ss(2)=s_flux_h_nw_2;
ss(3)=s_flux_v_nw_1;
ss(4)=s_flux_v_nw_2;
               
% ss(1)=sol_exactax_sigma(x2,y2,tt,nu,1,1);
% ss(2)=sol_exactax_sigma(x2,y2,tt,nu,2,1);
% ss(3)=sol_exactax_sigma(x2,y2,tt,nu,1,2);
% ss(4)=sol_exactax_sigma(x2,y2,tt,nu,2,2);

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
% sigma=sigma*(1/N);
return
end
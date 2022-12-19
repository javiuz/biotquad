function sigma=build_sigma_0_gral_grid(tt,nu)

global NN x y 

N=NN;

sigma=zeros(8*N*(N+1),1);

s1x=zeros(N+1,N+1,2);
s2x=s1x;
s1y=s1x;
s2y=s1x;

normal_n=s1x;
normal_s=normal_n;
normal_e=normal_n;
normal_w=normal_n;

% South-West corner node

% Comp. normales de las tensiones (1ª fila por un lado y 2ª fila por otro) 
% vienen dadas por normal_e y normal_n (celda noreste de los apuntes). 

vdim=4;
ss=zeros(vdim,1);

i=1;
j=1;

normal_n(i,j,1)=(y(i,j+1)-y(i,j))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
normal_n(i,j,2)=(x(i,j)-x(i,j+1))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
normal_e(i,j,1)=(y(i,j)-y(i+1,j))/sqrt((y(i,j)-y(i+1,j))^2+(x(i+1,j)-x(i,j))^2);
normal_e(i,j,2)=(x(i+1,j)-x(i,j))/sqrt((y(i,j)-y(i+1,j))^2+(x(i+1,j)-x(i,j))^2);

denom_ne=normal_e(i,j,1)*normal_n(i,j,2)-normal_n(i,j,1)*normal_e(i,j,2);
    alpha_ne=normal_n(i,j,2)/denom_ne;
    beta_ne=-normal_e(i,j,2)/denom_ne;
    gamma_ne=-normal_n(i,j,1)/denom_ne;
    delta_ne=normal_e(i,j,1)/denom_ne;

x1=x(i,j);  
y1=y(i,j);            
x2=x(i+1,j);
y2=y(i+1,j);
% x3=x(i+1,j+1);
% y3=y(i+1,j+1);
x4=x(i,j+1);
y4=y(i,j+1);

% We need to multiply to the module squared in order to obtain the
% fluxes multiplied by |e|.
mod_e1_sqrd=(x1-x2)^2+(y1-y2)^2;
mod_e2_sqrd=(x1-x4)^2+(y1-y4)^2;

s1y(i,j,1)=sol_exactax_sigma(x1,y1,tt,nu,1,2)*mod_e1_sqrd;
s2y(i,j,1)=sol_exactax_sigma(x1,y1,tt,nu,2,2)*mod_e1_sqrd;
s1x(i,j,2)=sol_exactax_sigma(x1,y1,tt,nu,1,1)*mod_e2_sqrd;
s2x(i,j,2)=sol_exactax_sigma(x1,y1,tt,nu,2,1)*mod_e2_sqrd;

    % Celda NE
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
    s_flux_h_ne_1=alpha_ne*s1y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                  beta_ne*s1x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
        % Comp. vertical del flujo 
    s_flux_v_ne_1=gamma_ne*s1y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                  delta_ne*s1x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
    % 2ª fila de la matriz de tensiones
        % Comp. horizontal del flujo         
    s_flux_h_ne_2=alpha_ne*s2y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                  beta_ne*s2x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
        % Comp. vertical del flujo 
    s_flux_v_ne_2=gamma_ne*s2y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                  delta_ne*s2x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
        
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
   
    x1=x(i-1,j);
    y1=y(i-1,j);
    x2=x(i,j);
    y2=y(i,j);
    x3=x(i,j+1);
    y3=y(i,j+1);
%     x4=x(i-1,j+1);
%     y4=y(i-1,j+1);
    x5=x(i+1,j);
    y5=y(i+1,j);
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
    
    mod_e1_sqrd=(x2-x5)^2+(y2-y5)^2;
    mod_e2_sqrd=(x2-x3)^2+(y2-y3)^2;
    mod_e3_sqrd=(x2-x1)^2+(y2-y1)^2;    
        
    s1y(i,j,1)=sol_exactax_sigma(x2,y2,tt,nu,1,2)*mod_e1_sqrd;
    s2y(i,j,1)=sol_exactax_sigma(x2,y2,tt,nu,2,2)*mod_e1_sqrd;
    s1x(i,j,2)=sol_exactax_sigma(x2,y2,tt,nu,1,1)*mod_e2_sqrd;
    s2x(i,j,2)=sol_exactax_sigma(x2,y2,tt,nu,2,1)*mod_e2_sqrd;
    s1y(i,j,2)=sol_exactax_sigma(x2,y2,tt,nu,1,2)*mod_e3_sqrd;
    s2y(i,j,2)=sol_exactax_sigma(x2,y2,tt,nu,2,2)*mod_e3_sqrd;
    
    % Celda NE
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_ne_1=alpha_ne*s1y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                   beta_ne*s1x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
            % Comp. vertical del flujo 
     s_flux_v_ne_1=gamma_ne*s1y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                   delta_ne*s1x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);   
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_ne_2=alpha_ne*s2y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                   beta_ne*s2x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
            % Comp. vertical del flujo 
     s_flux_v_ne_2=gamma_ne*s2y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                   delta_ne*s2x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);   

    % Celda NO
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_nw_1=alpha_nw*s1x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   beta_nw*s1y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_nw_1=gamma_nw*s1x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   delta_nw*s1y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);   
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_nw_2=alpha_nw*s2x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   beta_nw*s2y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_nw_2=gamma_nw*s2x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   delta_nw*s2y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);   
               
    ss(1)=s_flux_v_ne_1;
    ss(2)=s_flux_v_ne_2;
    ss(3)=(s_flux_h_ne_1 + s_flux_h_nw_1)/2;
    ss(4)=(s_flux_h_ne_2 + s_flux_h_nw_2)/2;
    ss(5)=s_flux_v_nw_1;
    ss(6)=s_flux_v_nw_2;
    
    indr=ld+1:ld+vdim;
    sigma(indr)=ss;
    ld=ld+vdim;
end

% South-East corner node (j=1)
i=N+1;
vdim=4;
ss=zeros(vdim,1);

normal_n(i,j,1)=(y(i,j+1)-y(i,j))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
normal_n(i,j,2)=(x(i,j)-x(i,j+1))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
normal_w(i,j,1)=(y(i-1,j)-y(i,j))/sqrt((y(i-1,j)-y(i,j))^2+(x(i,j)-x(i-1,j))^2);
normal_w(i,j,2)=(x(i,j)-x(i-1,j))/sqrt((y(i-1,j)-y(i,j))^2+(x(i,j)-x(i-1,j))^2);

denom_nw=normal_n(i,j,1)*normal_w(i,j,2)-normal_w(i,j,1)*normal_n(i,j,2);
    alpha_nw=normal_w(i,j,2)/denom_nw;
    beta_nw=-normal_n(i,j,2)/denom_nw;
    gamma_nw=-normal_w(i,j,1)/denom_nw;
    delta_nw=normal_n(i,j,1)/denom_nw;

x1=x(i-1,j);  
y1=y(i-1,j);            
x2=x(i,j);
y2=y(i,j);
x3=x(i,j+1);
y3=y(i,j+1);
% x4=x(i-1,j+1);
% y4=y(i-1,j+1);

mod_e1_sqrd=(x2-x3)^2+(y2-y3)^2;
mod_e2_sqrd=(x2-x1)^2+(y2-y1)^2;

s1x(i,j,2)=sol_exactax_sigma(x2,y2,tt,nu,1,1)*mod_e1_sqrd;
s2x(i,j,2)=sol_exactax_sigma(x2,y2,tt,nu,2,1)*mod_e1_sqrd;
s1y(i,j,2)=sol_exactax_sigma(x2,y2,tt,nu,1,2)*mod_e2_sqrd;
s2y(i,j,2)=sol_exactax_sigma(x2,y2,tt,nu,2,2)*mod_e2_sqrd;

    % Celda NO
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_nw_1=alpha_nw*s1x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   beta_nw*s1y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_nw_1=gamma_nw*s1x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   delta_nw*s1y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);   
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_nw_2=alpha_nw*s2x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   beta_nw*s2y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_nw_2=gamma_nw*s2x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   delta_nw*s2y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);  

ss(1)=s_flux_h_nw_1;
ss(2)=s_flux_h_nw_2;
ss(3)=s_flux_v_nw_1;
ss(4)=s_flux_v_nw_2;

indr=ld+1:ld+vdim;
sigma(indr)=ss;
ld=ld+vdim;

for j=2:N

    % West nodes
    i=1;
        
    vdim=6;
    ss=zeros(vdim,1);
    
    x1=x(i,j-1);
    y1=y(i,j-1);
%     x2=x(i+1,j-1);
%     y2=y(i+1,j-1);
    x3=x(i+1,j);
    y3=y(i+1,j);
    x4=x(i,j);
    y4=y(i,j);
%     x5=x(i+1,j+1);
%     y5=y(i+1,j+1);
    x6=x(i,j+1);
    y6=y(i,j+1);
    
    normal_s(i,j,1)=(y(i,j)-y(i,j-1))/sqrt((y(i,j)-y(i,j-1))^2+(x(i,j-1)-x(i,j))^2);
    normal_s(i,j,2)=(x(i,j-1)-x(i,j))/sqrt((y(i,j)-y(i,j-1))^2+(x(i,j-1)-x(i,j))^2);
    normal_e(i,j,1)=(y(i,j)-y(i+1,j))/sqrt((y(i,j)-y(i+1,j))^2+(x(i+1,j)-x(i,j))^2);
    normal_e(i,j,2)=(x(i+1,j)-x(i,j))/sqrt((y(i,j)-y(i+1,j))^2+(x(i+1,j)-x(i,j))^2);
    normal_n(i,j,1)=(y(i,j+1)-y(i,j))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
    normal_n(i,j,2)=(x(i,j)-x(i,j+1))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
    
    denom_se=normal_s(i,j,1)*normal_e(i,j,2)-normal_e(i,j,1)*normal_s(i,j,2);
        alpha_se=normal_e(i,j,2)/denom_se;
        beta_se=-normal_s(i,j,2)/denom_se;
        gamma_se=-normal_e(i,j,1)/denom_se;
        delta_se=normal_s(i,j,1)/denom_se; 

    denom_ne=normal_e(i,j,1)*normal_n(i,j,2)-normal_n(i,j,1)*normal_e(i,j,2);
        alpha_ne=normal_n(i,j,2)/denom_ne;
        beta_ne=-normal_e(i,j,2)/denom_ne;
        gamma_ne=-normal_n(i,j,1)/denom_ne;
        delta_ne=normal_e(i,j,1)/denom_ne; 
        
    mod_e1_sqrd=(x4-x1)^2+(y4-y1)^2;
    mod_e2_sqrd=(x4-x3)^2+(y4-y3)^2;
    mod_e3_sqrd=(x4-x6)^2+(y4-y6)^2;

    s1x(i,j,1)=sol_exactax_sigma(x4,y4,tt,nu,1,1)*mod_e1_sqrd;
    s2x(i,j,1)=sol_exactax_sigma(x4,y4,tt,nu,2,1)*mod_e1_sqrd;
    s1y(i,j,1)=sol_exactax_sigma(x4,y4,tt,nu,1,2)*mod_e2_sqrd;
    s2y(i,j,1)=sol_exactax_sigma(x4,y4,tt,nu,2,2)*mod_e2_sqrd;
    s1x(i,j,2)=sol_exactax_sigma(x4,y4,tt,nu,1,1)*mod_e3_sqrd;
    s2x(i,j,2)=sol_exactax_sigma(x4,y4,tt,nu,2,1)*mod_e3_sqrd;
    
    % Celda SE
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_se_1=alpha_se*s1x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
                   beta_se*s1y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_se_1=gamma_se*s1x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
                   delta_se*s1y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);  
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_se_2=alpha_se*s2x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
                   beta_se*s2y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_se_2=gamma_se*s2x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
                   delta_se*s2y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);    
    
    % Celda NE
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_ne_1=alpha_ne*s1y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                   beta_ne*s1x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
            % Comp. vertical del flujo 
     s_flux_v_ne_1=gamma_ne*s1y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                   delta_ne*s1x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);   
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_ne_2=alpha_ne*s2y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                   beta_ne*s2x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
            % Comp. vertical del flujo 
     s_flux_v_ne_2=gamma_ne*s2y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                   delta_ne*s2x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);   
               
    ss(1)=s_flux_h_se_1;
    ss(2)=s_flux_h_se_2;
    ss(3)=(s_flux_v_ne_1 + s_flux_v_se_1)/2;
    ss(4)=(s_flux_v_ne_2 + s_flux_v_se_2)/2;
    ss(5)=s_flux_h_ne_1;
    ss(6)=s_flux_h_ne_2;

    indr=ld+1:ld+vdim;
    sigma(indr)=ss;
    ld=ld+vdim;
    
    % Central nodes 
    vdim=8;
    ss=zeros(vdim,1);
         
    for i=2:N
        
%     x1=x(i-1,j-1);
%     y1=y(i-1,j-1);
    x2=x(i,j-1);
    y2=y(i,j-1);
    x3=x(i,j);
    y3=y(i,j);
    x4=x(i-1,j);
    y4=y(i-1,j);
%     x5=x(i+1,j-1);
%     y5=y(i+1,j-1);
    x6=x(i+1,j);
    y6=y(i+1,j);
%     x7=x(i+1,j+1);
%     y7=y(i+1,j+1);
    x8=x(i,j+1);
    y8=y(i,j+1);
%     x9=x(i-1,j+1);
%     y9=y(i-1,j+1);

    normal_s(i,j,1)=(y(i,j)-y(i,j-1))/sqrt((y(i,j)-y(i,j-1))^2+(x(i,j-1)-x(i,j))^2);
    normal_s(i,j,2)=(x(i,j-1)-x(i,j))/sqrt((y(i,j)-y(i,j-1))^2+(x(i,j-1)-x(i,j))^2);
    normal_e(i,j,1)=(y(i,j)-y(i+1,j))/sqrt((y(i,j)-y(i+1,j))^2+(x(i+1,j)-x(i,j))^2);
    normal_e(i,j,2)=(x(i+1,j)-x(i,j))/sqrt((y(i,j)-y(i+1,j))^2+(x(i+1,j)-x(i,j))^2);
    normal_n(i,j,1)=(y(i,j+1)-y(i,j))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
    normal_n(i,j,2)=(x(i,j)-x(i,j+1))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
    normal_w(i,j,1)=(y(i-1,j)-y(i,j))/sqrt((y(i-1,j)-y(i,j))^2+(x(i,j)-x(i-1,j))^2);
    normal_w(i,j,2)=(x(i,j)-x(i-1,j))/sqrt((y(i-1,j)-y(i,j))^2+(x(i,j)-x(i-1,j))^2);
    
    denom_se=normal_s(i,j,1)*normal_e(i,j,2)-normal_e(i,j,1)*normal_s(i,j,2);
        alpha_se=normal_e(i,j,2)/denom_se;
        beta_se=-normal_s(i,j,2)/denom_se;
        gamma_se=-normal_e(i,j,1)/denom_se;
        delta_se=normal_s(i,j,1)/denom_se;        
    
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
    
    denom_sw=normal_w(i,j,1)*normal_s(i,j,2)-normal_s(i,j,1)*normal_w(i,j,2);
        alpha_sw=normal_s(i,j,2)/denom_sw;
        beta_sw=-normal_w(i,j,2)/denom_sw;
        gamma_sw=-normal_s(i,j,1)/denom_sw;
        delta_sw=normal_w(i,j,1)/denom_sw;

    mod_e1_sqrd=(x3-x2)^2+(y3-y2)^2;
    mod_e2_sqrd=(x3-x6)^2+(y3-y6)^2;
    mod_e3_sqrd=(x3-x8)^2+(y3-y8)^2;
    mod_e4_sqrd=(x3-x4)^2+(y3-y4)^2;
        
    s1x(i,j,1)=sol_exactax_sigma(x3,y3,tt,nu,1,1)*mod_e1_sqrd;
    s2x(i,j,1)=sol_exactax_sigma(x3,y3,tt,nu,2,1)*mod_e1_sqrd;
    s1y(i,j,1)=sol_exactax_sigma(x3,y3,tt,nu,1,2)*mod_e2_sqrd;
    s2y(i,j,1)=sol_exactax_sigma(x3,y3,tt,nu,2,2)*mod_e2_sqrd;
    s1x(i,j,2)=sol_exactax_sigma(x3,y3,tt,nu,1,1)*mod_e3_sqrd;
    s2x(i,j,2)=sol_exactax_sigma(x3,y3,tt,nu,2,1)*mod_e3_sqrd;
    s1y(i,j,2)=sol_exactax_sigma(x3,y3,tt,nu,1,2)*mod_e4_sqrd;
    s2y(i,j,2)=sol_exactax_sigma(x3,y3,tt,nu,2,2)*mod_e4_sqrd;
    
    % Celda SE
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_se_1=alpha_se*s1x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
                   beta_se*s1y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_se_1=gamma_se*s1x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
                   delta_se*s1y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);  
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_se_2=alpha_se*s2x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
                   beta_se*s2y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_se_2=gamma_se*s2x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
                   delta_se*s2y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);    
    
    % Celda NE
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_ne_1=alpha_ne*s1y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                   beta_ne*s1x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
            % Comp. vertical del flujo 
     s_flux_v_ne_1=gamma_ne*s1y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                   delta_ne*s1x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);   
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_ne_2=alpha_ne*s2y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                   beta_ne*s2x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
            % Comp. vertical del flujo 
     s_flux_v_ne_2=gamma_ne*s2y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
                   delta_ne*s2x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
               
    % Celda NO
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_nw_1=alpha_nw*s1x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   beta_nw*s1y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_nw_1=gamma_nw*s1x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   delta_nw*s1y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);   
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_nw_2=alpha_nw*s2x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   beta_nw*s2y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_nw_2=gamma_nw*s2x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   delta_nw*s2y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);      
               
    % Celda SO
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_sw_1=alpha_sw*s1y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
                   beta_sw*s1x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);
            % Comp. vertical del flujo 
     s_flux_v_sw_1=gamma_sw*s1y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
                   delta_sw*s1x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);  
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_sw_2=alpha_sw*s2y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
                   beta_sw*s2x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);
            % Comp. vertical del flujo 
     s_flux_v_sw_2=gamma_sw*s2y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
                   delta_sw*s2x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);    
               
    ss(1)=(s_flux_h_se_1 + s_flux_h_sw_1)/2;
    ss(2)=(s_flux_h_se_2 + s_flux_h_sw_2)/2;
    ss(3)=(s_flux_v_se_1 + s_flux_v_ne_1)/2;
    ss(4)=(s_flux_v_se_2 + s_flux_v_ne_2)/2;
    ss(5)=(s_flux_h_ne_1 + s_flux_h_nw_1)/2;
    ss(6)=(s_flux_h_ne_2 + s_flux_h_nw_2)/2;
    ss(7)=(s_flux_v_sw_1 + s_flux_v_nw_1)/2;
    ss(8)=(s_flux_v_sw_2 + s_flux_v_nw_2)/2;  

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
    x2=x(i,j-1);
    y2=y(i,j-1);
    x3=x(i,j);
    y3=y(i,j);
    x4=x(i-1,j);
    y4=y(i-1,j);
    x5=x(i,j+1);
    y5=y(i,j+1);
%     x6=x(i-1,j+1);
%     y6=y(i-1,j+1);   

    normal_s(i,j,1)=(y(i,j)-y(i,j-1))/sqrt((y(i,j)-y(i,j-1))^2+(x(i,j-1)-x(i,j))^2);
    normal_s(i,j,2)=(x(i,j-1)-x(i,j))/sqrt((y(i,j)-y(i,j-1))^2+(x(i,j-1)-x(i,j))^2);
    normal_n(i,j,1)=(y(i,j+1)-y(i,j))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
    normal_n(i,j,2)=(x(i,j)-x(i,j+1))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
    normal_w(i,j,1)=(y(i-1,j)-y(i,j))/sqrt((y(i-1,j)-y(i,j))^2+(x(i,j)-x(i-1,j))^2);
    normal_w(i,j,2)=(x(i,j)-x(i-1,j))/sqrt((y(i-1,j)-y(i,j))^2+(x(i,j)-x(i-1,j))^2);       
    
    denom_nw=normal_n(i,j,1)*normal_w(i,j,2)-normal_w(i,j,1)*normal_n(i,j,2);
        alpha_nw=normal_w(i,j,2)/denom_nw;
        beta_nw=-normal_n(i,j,2)/denom_nw;
        gamma_nw=-normal_w(i,j,1)/denom_nw;
        delta_nw=normal_n(i,j,1)/denom_nw;
    
    denom_sw=normal_w(i,j,1)*normal_s(i,j,2)-normal_s(i,j,1)*normal_w(i,j,2);
        alpha_sw=normal_s(i,j,2)/denom_sw;
        beta_sw=-normal_w(i,j,2)/denom_sw;
        gamma_sw=-normal_s(i,j,1)/denom_sw;
        delta_sw=normal_w(i,j,1)/denom_sw;
    
    mod_e1_sqrd=(x3-x2)^2+(y3-y2)^2;
    mod_e2_sqrd=(x3-x5)^2+(y3-y5)^2;
    mod_e3_sqrd=(x3-x4)^2+(y3-y4)^2;
        
    s1x(i,j,1)=sol_exactax_sigma(x3,y3,tt,nu,1,1)*mod_e1_sqrd;
    s2x(i,j,1)=sol_exactax_sigma(x3,y3,tt,nu,2,1)*mod_e1_sqrd;
    s1x(i,j,2)=sol_exactax_sigma(x3,y3,tt,nu,1,1)*mod_e2_sqrd;
    s2x(i,j,2)=sol_exactax_sigma(x3,y3,tt,nu,2,1)*mod_e2_sqrd;
    s1y(i,j,2)=sol_exactax_sigma(x3,y3,tt,nu,1,2)*mod_e3_sqrd;
    s2y(i,j,2)=sol_exactax_sigma(x3,y3,tt,nu,2,2)*mod_e3_sqrd;
    
    % Celda NO
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_nw_1=alpha_nw*s1x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   beta_nw*s1y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_nw_1=gamma_nw*s1x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   delta_nw*s1y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);   
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_nw_2=alpha_nw*s2x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   beta_nw*s2y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_nw_2=gamma_nw*s2x(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
                   delta_nw*s2y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);      
               
    % Celda SO
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_sw_1=alpha_sw*s1y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
                   beta_sw*s1x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);
            % Comp. vertical del flujo 
     s_flux_v_sw_1=gamma_sw*s1y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
                   delta_sw*s1x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);  
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_sw_2=alpha_sw*s2y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
                   beta_sw*s2x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);
            % Comp. vertical del flujo 
     s_flux_v_sw_2=gamma_sw*s2y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
                   delta_sw*s2x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);  
               
    ss(1)=s_flux_h_sw_1;
    ss(2)=s_flux_h_sw_2;
    ss(3)=s_flux_h_nw_1;
    ss(4)=s_flux_h_nw_2;
    ss(5)=(s_flux_v_sw_1 + s_flux_v_nw_1)/2;  
    ss(6)=(s_flux_v_sw_2 + s_flux_v_nw_2)/2;            

    indr=ld+1:ld+vdim;
    sigma(indr)=ss;
    ld=ld+vdim;
end

% North-West corner node
i=1;
vdim=4;
ss=zeros(vdim,1);

x1=x(i,j-1);
y1=y(i,j-1);           
% x2=x(i+1,j-1);
% y2=y(i+1,j-1);
x3=x(i+1,j);
y3=y(i+1,j);
x4=x(i,j);
y4=y(i,j);

denom_se=normal_s(i,j,1)*normal_e(i,j,2)-normal_e(i,j,1)*normal_s(i,j,2);
    alpha_se=normal_e(i,j,2)/denom_se;
    beta_se=-normal_s(i,j,2)/denom_se;
    gamma_se=-normal_e(i,j,1)/denom_se;
    delta_se=normal_s(i,j,1)/denom_se;   
    
mod_e1_sqrd=(x4-x1)^2+(y4-y1)^2;
mod_e2_sqrd=(x4-x3)^2+(y4-y3)^2;
    
s1x(i,j,1)=sol_exactax_sigma(x4,y4,tt,nu,1,1)*mod_e1_sqrd;
s2x(i,j,1)=sol_exactax_sigma(x4,y4,tt,nu,2,1)*mod_e1_sqrd;
s1y(i,j,1)=sol_exactax_sigma(x4,y4,tt,nu,1,2)*mod_e2_sqrd;
s2y(i,j,1)=sol_exactax_sigma(x4,y4,tt,nu,2,2)*mod_e2_sqrd;
    
    % Celda SE
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_se_1=alpha_se*s1x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
                   beta_se*s1y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_se_1=gamma_se*s1x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
                   delta_se*s1y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);  
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_se_2=alpha_se*s2x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
                   beta_se*s2y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_se_2=gamma_se*s2x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
                   delta_se*s2y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);    
               
ss(1)=s_flux_h_se_1;
ss(2)=s_flux_h_se_2;
ss(3)=s_flux_v_se_1;
ss(4)=s_flux_v_se_2;

indr=ld+1:ld+vdim;
sigma(indr)=ss;
ld=ld+vdim;

% North nodes (j=N+1)
vdim=6;
ss=zeros(vdim,1);

for i=2:N
    
%     x1=x(i-1,j-1);
%     y1=y(i-1,j-1);
    x2=x(i,j-1);
    y2=y(i,j-1);
    x3=x(i,j);
    y3=y(i,j);
    x4=x(i-1,j);
    y4=y(i-1,j);
%     x5=x(i+1,j-1);
%     y5=y(i+1,j-1);
    x6=x(i+1,j);
    y6=y(i+1,j);   

    normal_s(i,j,1)=(y(i,j)-y(i,j-1))/sqrt((y(i,j)-y(i,j-1))^2+(x(i,j-1)-x(i,j))^2);
    normal_s(i,j,2)=(x(i,j-1)-x(i,j))/sqrt((y(i,j)-y(i,j-1))^2+(x(i,j-1)-x(i,j))^2);
    normal_e(i,j,1)=(y(i,j)-y(i+1,j))/sqrt((y(i,j)-y(i+1,j))^2+(x(i+1,j)-x(i,j))^2);
    normal_e(i,j,2)=(x(i+1,j)-x(i,j))/sqrt((y(i,j)-y(i+1,j))^2+(x(i+1,j)-x(i,j))^2);
    normal_w(i,j,1)=(y(i-1,j)-y(i,j))/sqrt((y(i-1,j)-y(i,j))^2+(x(i,j)-x(i-1,j))^2);
    normal_w(i,j,2)=(x(i,j)-x(i-1,j))/sqrt((y(i-1,j)-y(i,j))^2+(x(i,j)-x(i-1,j))^2);
    
    denom_se=normal_s(i,j,1)*normal_e(i,j,2)-normal_e(i,j,1)*normal_s(i,j,2);
        alpha_se=normal_e(i,j,2)/denom_se;
        beta_se=-normal_s(i,j,2)/denom_se;
        gamma_se=-normal_e(i,j,1)/denom_se;
        delta_se=normal_s(i,j,1)/denom_se;        
    
    denom_sw=normal_w(i,j,1)*normal_s(i,j,2)-normal_s(i,j,1)*normal_w(i,j,2);
        alpha_sw=normal_s(i,j,2)/denom_sw;
        beta_sw=-normal_w(i,j,2)/denom_sw;
        gamma_sw=-normal_s(i,j,1)/denom_sw;
        delta_sw=normal_w(i,j,1)/denom_sw;

    mod_e1_sqrd=(x3-x2)^2+(y3-y2)^2;
    mod_e2_sqrd=(x3-x6)^2+(y3-y6)^2;
    mod_e3_sqrd=(x3-x4)^2+(y3-y4)^2;
        
    s1x(i,j,1)=sol_exactax_sigma(x3,y3,tt,nu,1,1)*mod_e1_sqrd;
    s2x(i,j,1)=sol_exactax_sigma(x3,y3,tt,nu,2,1)*mod_e1_sqrd;
    s1y(i,j,1)=sol_exactax_sigma(x3,y3,tt,nu,1,2)*mod_e2_sqrd;
    s2y(i,j,1)=sol_exactax_sigma(x3,y3,tt,nu,2,2)*mod_e2_sqrd;
    s1y(i,j,2)=sol_exactax_sigma(x3,y3,tt,nu,1,2)*mod_e3_sqrd;
    s2y(i,j,2)=sol_exactax_sigma(x3,y3,tt,nu,2,2)*mod_e3_sqrd;
    
    % Celda SE
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_se_1=alpha_se*s1x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
                   beta_se*s1y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_se_1=gamma_se*s1x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
                   delta_se*s1y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);  
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_se_2=alpha_se*s2x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
                   beta_se*s2y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);
            % Comp. vertical del flujo 
     s_flux_v_se_2=gamma_se*s2x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
                   delta_se*s2y(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);       
               
    % Celda SO
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_sw_1=alpha_sw*s1y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
                   beta_sw*s1x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);
            % Comp. vertical del flujo 
     s_flux_v_sw_1=gamma_sw*s1y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
                   delta_sw*s1x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);  
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_sw_2=alpha_sw*s2y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
                   beta_sw*s2x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);
            % Comp. vertical del flujo 
     s_flux_v_sw_2=gamma_sw*s2y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
                   delta_sw*s2x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);  
               
    ss(1)=(s_flux_h_se_1 + s_flux_h_sw_1)/2;
    ss(2)=(s_flux_h_se_2 + s_flux_h_sw_2)/2;
    ss(3)=s_flux_v_se_1;
    ss(4)=s_flux_v_se_2;
    ss(5)=s_flux_v_sw_1;
    ss(6)=s_flux_v_sw_2;

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
x2=x(i,j-1);
y2=y(i,j-1);
x3=x(i,j);
y3=y(i,j);
x4=x(i-1,j);
y4=y(i-1,j);

denom_sw=normal_w(i,j,1)*normal_s(i,j,2)-normal_s(i,j,1)*normal_w(i,j,2);
    alpha_sw=normal_s(i,j,2)/denom_sw;
    beta_sw=-normal_w(i,j,2)/denom_sw;
    gamma_sw=-normal_s(i,j,1)/denom_sw;
    delta_sw=normal_w(i,j,1)/denom_sw;
    
mod_e1_sqrd=(x3-x2)^2+(y3-y2)^2;
mod_e2_sqrd=(x3-x4)^2+(y3-y4)^2;
    
s1x(i,j,1)=sol_exactax_sigma(x3,y3,tt,nu,1,1)*mod_e1_sqrd;
s2x(i,j,1)=sol_exactax_sigma(x3,y3,tt,nu,2,1)*mod_e1_sqrd;
s1y(i,j,2)=sol_exactax_sigma(x3,y3,tt,nu,1,2)*mod_e2_sqrd;
s2y(i,j,2)=sol_exactax_sigma(x3,y3,tt,nu,2,2)*mod_e2_sqrd;
    
    % Celda SO
        % 1ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_sw_1=alpha_sw*s1y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
                   beta_sw*s1x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);
            % Comp. vertical del flujo 
     s_flux_v_sw_1=gamma_sw*s1y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
                   delta_sw*s1x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);  
        % 2ª fila de la matriz de tensiones
            % Comp. horizontal del flujo 
     s_flux_h_sw_2=alpha_sw*s2y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
                   beta_sw*s2x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);
            % Comp. vertical del flujo 
     s_flux_v_sw_2=gamma_sw*s2y(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
                   delta_sw*s2x(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);      

ss(1)=s_flux_h_sw_1;
ss(2)=s_flux_h_sw_2;
ss(3)=s_flux_v_sw_1;
ss(4)=s_flux_v_sw_2;      

indr=ld+1:ld+vdim;
sigma(indr)=ss;
return
end
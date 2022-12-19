function sigma=build_sigma_0_gral_grid_3(s1x,s1y,s2x,s2y)

global NN x y 

N=NN;

sigma=zeros(8*N*(N+1),1);

normal_n=zeros(N+1,N+1,2);
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

indr=1:vdim;
sigma(indr)=ss;
ld=vdim;

% South nodes (j=1)
vdim=6;
ss=zeros(vdim,1);

for i=2:N

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

denom_se=normal_s(i,j,1)*normal_e(i,j,2)-normal_e(i,j,1)*normal_s(i,j,2);
    alpha_se=normal_e(i,j,2)/denom_se;
    beta_se=-normal_s(i,j,2)/denom_se;
    gamma_se=-normal_e(i,j,1)/denom_se;
    delta_se=normal_s(i,j,1)/denom_se;    
    
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

denom_sw=normal_w(i,j,1)*normal_s(i,j,2)-normal_s(i,j,1)*normal_w(i,j,2);
    alpha_sw=normal_s(i,j,2)/denom_sw;
    beta_sw=-normal_w(i,j,2)/denom_sw;
    gamma_sw=-normal_s(i,j,1)/denom_sw;
    delta_sw=normal_w(i,j,1)/denom_sw;
    
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
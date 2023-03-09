function MS_MF_MFE(N)

global NN x y lambda mu alpha c0 perm
 
NN=N;       % Dimensión del problema discreto.

% Generación de la malla 
mesh=3;                 
[x,y]=init_mesh(mesh);  % coordenadas de los vértices de la malla.

% Parámetros ecuación Biot
% alpha=1;
alpha=0;
lambda=0;
mu=1;
% lambda=123;
% mu=79.3;

% Parámetros ecuación flujo
perm=1;

% Storativity coefficient
c0=0;
% c0=1e-05;

% Hydraulic conductivity: is inside the function 'kinv.m'

% Time step
t=0;
delta_t=1e-04;
Tf=1e-03;

%  Matriz de fuerza: stress-stress
A_ss=build_matrix_sig_sig;

% Matriz de fuerza (transpuesta): stress-displacement
A_su_T=build_matrix_sig_u;

% Matriz de fuerza: A_sigma_u
A_su=(A_su_T)';

% Matriz de fuerza (transpuesta): stress-rotation
A_sg_T=build_matrix_sig_gamma;

% Matriz de fuerza: A_sigma_gamma
A_sg=(A_sg_T)';

% Matriz de fuerza (transpuesta): stress-pressure
A_sp_T=build_matrix_p_tau;
% A_sp_T=zeros(8*N*(N+1),N*N);

% Matriz de fuerza: A_sigma_p
A_sp=(A_sp_T)';
% A_sp=zeros(N*N,8*N*(N+1));

% Matriz de fuerza: pressure-pressure
A_pp=(c0+alpha^2/(lambda+mu))*(1/(N^2))*eye(N*N);
% A_pp=0*eye(N*N);

% Matriz de fuerza: velocity-velocity
A_zz=build_matrix_z_z;

% Matriz de fuerza (transpuesta): velocity-pressure
A_zp_T=build_matrix_z_p;

% Matriz de fuerza: A_zp
A_zp=(A_zp_T)'*delta_t;

M_big=[A_ss A_su_T A_sg_T zeros(8*N*(N+1),4*N*(N+1)) A_sp_T;...
       -A_su zeros(2*N*N,2*N*N+(N+1)*(N+1)+4*N*(N+1)+N*N);...
       -A_sg zeros((N+1)*(N+1),2*N*N+(N+1)*(N+1)+4*N*(N+1)+N*N);...
       zeros(4*N*(N+1),8*N*(N+1)+2*N*N+(N+1)*(N+1)) A_zz A_zp_T;...
       A_sp zeros(N*N,2*N*N+(N+1)*(N+1)) -A_zp A_pp];
   
% eigs(M_big)
% pause
   
% M_big_st=[A_ss A_su_T A_sg_T zeros(8*N*(N+1),4*N*(N+1)) A_sp_T;...
%        -A_su zeros(2*N*N,2*N*N+(N+1)*(N+1)+4*N*(N+1)+N*N);...
%        -A_sg zeros((N+1)*(N+1),2*N*N+(N+1)*(N+1)+4*N*(N+1)+N*N);...
%        zeros(4*N*(N+1),8*N*(N+1)+2*N*N+(N+1)*(N+1)) A_zz A_zp_T;...
%        zeros(N*N,8*N*(N+1)+2*N*N+(N+1)*(N+1)) -A_zp zeros(N*N,N*N)];

   % Initial solution of the variables at t=0 
% u=zeros(2*N*N,1);
p=zeros(N*N,1);

for j=1:N
    for i=1:N
        ind2u=(i+(j-1)*N)*2;
%         ind1u=ind2u-1;
        ind1p=ind2u/2;
        
        xx=(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1))/4;
        yy=(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1))/4;
        
%         u(ind1u)=sol_exactax(xx,yy,t,1);
%         u(ind2u)=sol_exactax(xx,yy,t,2);
        p(ind1p)=sol_exactax(xx,yy,t,3);
    end
end

nu=0;
% sigma=build_sigma_0_cartesian(t,nu);
sigma=build_sigma_0_gral_grid(t,nu);

%% Initialize errors
erroru_L2_inf=0;
erroru_3_inf=erroru_L2_inf;
errorg_L2_inf=erroru_L2_inf;
error_sigma_sigmah_inf=erroru_L2_inf;
errorp_L2_inf=erroru_L2_inf;
errorp_3_inf=erroru_L2_inf;
error_z_zh_inf=erroru_L2_inf;

% indep_old=zeros(12*N*(N+1)+3*N*N+(N+1)*(N+1),1);
indep1old=zeros(8*N*(N+1),1);
indep2old=zeros(2*N*N,1);
indep4old=zeros(4*N*(N+1),1);
indep5old=zeros(N*N,1);

%% We solve Biot's system: Time loop

while t < Tf

indep2=build_indep_f(t+delta_t);
% indep2=0*indep2;

indep5=delta_t*build_indep_q(t+delta_t) + A_sp*sigma + A_pp*p;

[indep1,indep4]=dir_bc_Pg(t+delta_t);

indep=[indep1;indep2;zeros((N+1)*(N+1),1);indep4;indep5];
   
Sol=M_big\indep;

% disp(max(abs(indep-indep_old)))

% disp(max(abs(indep1-indep1old)))
% disp(max(abs(indep2-indep2old)))
% disp(max(abs(indep4-indep4old)))
% disp(max(abs(indep5-indep5old)))
% pause
% indep1old=indep1;
% indep2old=indep2;
% indep4old=indep4;
% indep5old=indep5;

% indep_old=indep;

sigma=Sol(1:8*N*(N+1));
u=Sol(8*N*(N+1)+1:8*N*(N+1)+2*N*N);
gamma=Sol(8*N*(N+1)+2*N*N+1:8*N*(N+1)+2*N*N+(N+1)*(N+1));
z=Sol(8*N*(N+1)+2*N*N+(N+1)*(N+1)+1:12*N*(N+1)+2*N*N+(N+1)*(N+1));
p=Sol(12*N*(N+1)+2*N*N+(N+1)*(N+1)+1:12*N*(N+1)+3*N*N+(N+1)*(N+1));

% Update t
t=t+delta_t;

    % We reorder the variables to compute the errors and some contourplots
gamma_n=reshape(gamma,N+1,N+1);
sigma_n=build_sigma_n_2(sigma);
vel_n=build_vel_n(z);

    % L2 norms of the different variables
% t
        % Displacement
[erroru_L2,erroru_3]=compute_error_displacements(u,t);
        % Pressure
[errorp_L2,errorp_3]=compute_error_pressures(p,t);
        % Rotation
errorg_L2=compute_error_rotations(gamma_n,t);
        % Stress
error_sigma_sigmah=compute_error_stress(sigma_n,nu,t);
        % Velocity
error_z_zh=compute_error_velocities(vel_n,t);

    % Infinity norm of the L2 errors of the variables
erroru_L2_inf=max(erroru_L2_inf,erroru_L2);
erroru_3_inf=max(erroru_3_inf,erroru_3);
errorp_L2_inf=max(errorp_L2_inf,errorp_L2);
errorp_3_inf=max(errorp_3_inf,errorp_3);
errorg_L2_inf=max(errorg_L2_inf,errorg_L2);
error_sigma_sigmah_inf=max(error_sigma_sigmah_inf,error_sigma_sigmah);
error_z_zh_inf=max(error_z_zh_inf,error_z_zh);
end

% We display the infinity norms of the errors of the variables at final time step
erroru_L2_inf
erroru_3_inf
errorp_L2_inf
errorp_3_inf
errorg_L2_inf
error_sigma_sigmah_inf
error_z_zh_inf

%% Plots of the Magnitude/Vectors of the variables

% % We associate the values of sigma_n to the vertices of the mesh.
% [s1x,s1y,s2x,s2y]=reorder_sigma_n(sigma_n);
% [stress_flux1]=compute_stress_fluxes_2(s1x,s1y);
% [stress_flux2]=compute_stress_fluxes_2(s2x,s2y);
% 
% [vel_flux]=compute_vel_fluxes(zx,zy);
% 
% Uc=sparse(N,N);
% Gc=Uc;
% Sc1=Uc;
% Sc2=Uc;
% Zc=Uc;
% Pc=reshape(p,N,N);
% Xc=Uc;
% Yc=Uc;
% U_ex=Uc;
% G_ex=Uc;
% S_ex1=Uc;
% S_ex2=Uc;
% Z_ex=Uc;
% P_ex=Uc;
% 
% for j=1:N
%     for i=1:N
%         ind2u=(i+(j-1)*N)*2;
%         ind1u=ind2u-1;
%         ind1p=ind2u/2;
%         xx=(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1))/4;
%         yy=(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1))/4;
%         Xc(i,j)=xx;
%         Yc(i,j)=yy;
%         U_ex(i,j)=sqrt(sol_exactax(xx,yy,t,1)^(2)+sol_exactax(xx,yy,t,2)^(2));
%         P_ex(i,j)=sol_exactax(xx,yy,t,3);
%         Z_ex(i,j)=sqrt(sol_exactax(xx,yy,t,4)^(2)+sol_exactax(xx,yy,t,5)^(2));
%         G_ex(i,j)=sol_exactax(xx,yy,t,6);
%         S_ex1(i,j)=sqrt(sol_exactax_sigma(xx,yy,t,nu,1,1)^(2)+sol_exactax_sigma(xx,yy,t,nu,1,2)^(2));
%         S_ex2(i,j)=sqrt(sol_exactax_sigma(xx,yy,t,nu,2,1)^(2)+sol_exactax_sigma(xx,yy,t,nu,2,2)^(2));
%         Uc(i,j)=sqrt((u(ind1u))^(2)+(u(ind2u))^(2));
%         Gc(i,j)=(gamma_n(i,j)+gamma_n(i+1,j)+gamma_n(i+1,j+1)+gamma_n(i,j+1))/4;
%         Sc1(i,j)=(sqrt(stress_flux1(i,j,1)^2+stress_flux1(i,j,2)^2)+...
%                   sqrt(stress_flux1(i+1,j,1)^2+stress_flux1(i+1,j,2)^2)+...
%                   sqrt(stress_flux1(i+1,j+1,1)^2+stress_flux1(i+1,j+1,2)^2)+...
%                   sqrt(stress_flux1(i,j+1,1)^2+stress_flux1(i,j+1,2)^2))/4;
%         Sc2(i,j)=(sqrt(stress_flux2(i,j,1)^2+stress_flux2(i,j,2)^2)+...
%                   sqrt(stress_flux2(i+1,j,1)^2+stress_flux2(i+1,j,2)^2)+...
%                   sqrt(stress_flux2(i+1,j+1,1)^2+stress_flux2(i+1,j+1,2)^2)+...
%                   sqrt(stress_flux2(i,j+1,1)^2+stress_flux2(i,j+1,2)^2))/4;
%         Zc(i,j)=(sqrt(vel_flux(i,j,1)^2+vel_flux(i,j,2)^2)+...
%                  sqrt(vel_flux(i+1,j,1)^2+vel_flux(i+1,j,2)^2)+...
%                  sqrt(vel_flux(i+1,j+1,1)^2+vel_flux(i+1,j+1,2)^2)+...
%                  sqrt(vel_flux(i,j+1,1)^2+vel_flux(i,j+1,2)^2))/4;
%     end
% end
% 
% %       p -> contourplot
% 
% figure
% contourf(Xc,Yc,P_ex); % contourlines for the pressure (exact solution)
% axis equal		
% % shading flat;
% grid off		
% set(gca, 'YTick', []);	
% set(gca, 'XTick', []);
% colorbar
% title('Darcy Pressure:sol exact')
% 
% figure
% contourf(Xc,Yc,Pc); % contourlines for the pressure (computed solution)
% axis equal		
% % shading flat;
% grid off		
% set(gca, 'YTick', []);	
% set(gca, 'XTick', []);
% colorbar
% title('Darcy Pressure:computed sol')
% 
% %       u -> contourplot
% 
% figure
% contourf(Xc,Yc,U_ex); % contourlines for the displacement (exact solution)
% axis equal		
% % shading flat;
% grid off		
% set(gca, 'YTick', []);	
% set(gca, 'XTick', []);
% colorbar
% title('Displacement (Euclidean norm):sol exact')
% 
% figure
% contourf(Xc,Yc,Uc); % contourlines for the displacement (computed solution)
% axis equal		
% % shading flat;
% grid off		
% set(gca, 'YTick', []);	
% set(gca, 'XTick', []);
% colorbar
% title('Displacement (Euclidean norm):computed sol')
% 
% %       gamma -> contourplot
% 
% figure
% contourf(Xc,Yc,G_ex); % contourlines for the rotation (exact solution)
% axis equal		
% % shading flat;
% grid off		
% set(gca, 'YTick', []);	
% set(gca, 'XTick', []);
% colorbar
% title('Rotation: sol exact')
% 
% figure
% contourf(Xc,Yc,Gc); % contourlines for the rotation (computed solution)
% axis equal		
% % shading flat;
% grid off		
% set(gca, 'YTick', []);	
% set(gca, 'XTick', []);
% colorbar
% title('Rotation: computed sol')
% 
% %       z -> contourplot
% 
% figure
% contourf(Xc,Yc,Z_ex); % contourlines for the velocity (exact solution)
% axis equal		
% % shading flat;
% grid off		
% set(gca, 'YTick', []);	
% set(gca, 'XTick', []);
% colorbar
% title('Darcy velocity (Euclidean norm):sol exact')
% 
% figure
% contourf(Xc,Yc,Zc); % contourlines for the velocity (computed solution)
% axis equal		
% % shading flat;
% grid off		
% set(gca, 'YTick', []);	
% set(gca, 'XTick', []);
% colorbar
% title('Darcy velocity (Euclidean norm):computed sol')
% 
% %       sigma (1st row)-> contourplot
% 
% figure
% contourf(Xc,Yc,S_ex1); % contourlines for the 1st row of the stress (exact solution)
% axis equal		
% % shading flat;
% grid off		
% set(gca, 'YTick', []);	
% set(gca, 'XTick', []);
% colorbar
% title('Stress (1st row):sol exact')   
%     
% figure
% contourf(Xc,Yc,Sc1); % contourlines for the 1st row of the stress (computed solution)
% axis equal		
% % shading flat;
% grid off		
% set(gca, 'YTick', []);	
% set(gca, 'XTick', []);
% colorbar
% title('Stress (1st row):computed sol')
% 
% figure
% contourf(Xc,Yc,S_ex2); % contourlines for the 2nd row of the stress (exact solution)
% axis equal		
% % shading flat;
% grid off		
% set(gca, 'YTick', []);	
% set(gca, 'XTick', []);
% colorbar
% title('Stress (2nd row):sol exact')    
%     
% figure
% contourf(Xc,Yc,Sc2); % contourlines for the 2nd row of the stress (computed sol)
% axis equal		
% % shading flat;
% grid off		
% set(gca, 'YTick', []);	
% set(gca, 'XTick', []);
% colorbar
% title('Stress (2nd row):computed sol')

end
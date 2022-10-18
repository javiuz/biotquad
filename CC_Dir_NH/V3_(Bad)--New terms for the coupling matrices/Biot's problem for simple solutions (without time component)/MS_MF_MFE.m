function MS_MF_MFE(N)

global NN x y alpha lambda mu c0 perm

NN=N;       % Dimensi�n del problema discreto.

% Generaci�n de la malla: �OJO! Para la malla n� 3 tenemos que introducir 
% el valor del n� de refinamientos de manera manual.
mesh=0;                 
[x,y]=init_mesh(mesh);  % coordenadas de los v�rtices de la malla.
tic

% Par�metros ecuaci�n Biot
alpha=1;
E=1;
nu=0.2;
lambda=(E*nu)/((1+nu)*(1-2*nu));
mu=E/(2*(1+nu));

% Storativity coefficient
c0=1e-05;

% Hydraulic conductivity: is inside the function 'kinv.m'
% K=perm*[1 0;0 1]; with perm=1, 1e-03, 1e-06, 1e-09, 1e-12
perm=1;

% % initial time
% t=0;
% % Time step
% delta_t=1e-04;
% % Final time
% Tf=1e-03;

%  Matrices del sistema de Biot: A11, A12, A21 y A22
    % Asp y App las utilizaremos despu�s
[A11,A12,A22]=build_matrices_Biot;
A21=zeros(N*N,2*N*N);

% Matriz del sistema reducido de Biot
Biot_matrix=[A11 A12;A21 A22];

% Source terms of the MFMFE-MSMFE discretization at t+delta_t
f_indep=build_indep_f;  % Source term f
% q_indep=zeros(N*N,1);   % Source term q

    % Non-homogeneous Dir. B.C. for u and p
[gDu,gDp]=dir_bc_Pg;
f_hat= -f_indep + gDu; 
q_hat= gDp;

% Right-hand side of the Biot system
indep_term=[f_hat;q_hat];

% We solve Biot's system
sol_vec=Biot_matrix\indep_term; 
u=sol_vec(1:2*N*N);
p=sol_vec(2*N*N+1:2*N*N+N*N);

% We compute the rest of the variables:
    % rotation 
gamma=compute_gamma(u,p);  % Computed solution for the rotation term
    % stress
[sigma,~,~,~,~]=compute_tensors(u,p,gamma);
    % velocity
[z,zx,zy]=compute_fluxes(p); 

% We reorder the variables to compute the errors and some contourplots
gamma_n=reshape(gamma,N+1,N+1);
sigma_n=build_sigma_n_2(sigma);
vel_n=build_vel_n(z);

% L2 norms of the different variables
    % Displacement
[erroru_L2,erroru_3]=compute_error_displacements(u)
    % Pressure
[errorp_L2,errorp_3]=compute_error_pressures(p)
    % Rotation
errorg_L2=compute_error_rotations(gamma_n)
    % Stress
error_sigma_sigmah=compute_error_stress(sigma_n,nu)
    % Velocity
error_z_zh=compute_error_velocities(vel_n)

% Components of the displacement vector: u1 and u2
% u1=u(1:2:2*N*N);
% u2=u(2:2:2*N*N);

%% Plots of the Magnitude/Vectors of the variables

% We associate the values of sigma_n to the vertices of the mesh.
[s1x,s1y,s2x,s2y]=reorder_sigma_n(sigma_n);
[stress_flux1]=compute_stress_fluxes_2(s1x,s1y);
[stress_flux2]=compute_stress_fluxes_2(s2x,s2y);

[vel_flux]=compute_vel_fluxes(zx,zy);

Uc=sparse(N,N);
Gc=Uc;
Sc1=Uc;
Sc2=Uc;
Zc=Uc;
Pc=reshape(p,N,N);
Xc=Uc;
Yc=Uc;
U_ex=Uc;
G_ex=Uc;
S_ex1=Uc;
S_ex2=Uc;
Z_ex=Uc;
P_ex=Uc;

for j=1:N
    for i=1:N
        ind2u=(i+(j-1)*N)*2;
        ind1u=ind2u-1;
%         ind1p=ind2u/2;
        xx=(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1))/4;
        yy=(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1))/4;
        Xc(i,j)=xx;
        Yc(i,j)=yy;
        U_ex(i,j)=sqrt(sol_exactax(xx,yy,1)^(2)+sol_exactax(xx,yy,2)^(2));
        P_ex(i,j)=sol_exactax(xx,yy,3);
        Z_ex(i,j)=sqrt(sol_exactax(xx,yy,4)^(2)+sol_exactax(xx,yy,5)^(2));
        G_ex(i,j)=sol_exactax(xx,yy,6);
        S_ex1(i,j)=sqrt(sol_exactax_sigma(xx,yy,nu,1,1)^(2)+sol_exactax_sigma(xx,yy,nu,1,2)^(2));
        S_ex2(i,j)=sqrt(sol_exactax_sigma(xx,yy,nu,2,1)^(2)+sol_exactax_sigma(xx,yy,nu,2,2)^(2));
        Uc(i,j)=sqrt((u(ind1u))^(2)+(u(ind2u))^(2));
        Gc(i,j)=(gamma_n(i,j)+gamma_n(i+1,j)+gamma_n(i+1,j+1)+gamma_n(i,j+1))/4;
        Sc1(i,j)=(sqrt(stress_flux1(i,j,1)^2+stress_flux1(i,j,2)^2)+...
                  sqrt(stress_flux1(i+1,j,1)^2+stress_flux1(i+1,j,2)^2)+...
                  sqrt(stress_flux1(i+1,j+1,1)^2+stress_flux1(i+1,j+1,2)^2)+...
                  sqrt(stress_flux1(i,j+1,1)^2+stress_flux1(i,j+1,2)^2))/4;
        Sc2(i,j)=(sqrt(stress_flux2(i,j,1)^2+stress_flux2(i,j,2)^2)+...
                  sqrt(stress_flux2(i+1,j,1)^2+stress_flux2(i+1,j,2)^2)+...
                  sqrt(stress_flux2(i+1,j+1,1)^2+stress_flux2(i+1,j+1,2)^2)+...
                  sqrt(stress_flux2(i,j+1,1)^2+stress_flux2(i,j+1,2)^2))/4;
        Zc(i,j)=(sqrt(vel_flux(i,j,1)^2+vel_flux(i,j,2)^2)+...
                 sqrt(vel_flux(i+1,j,1)^2+vel_flux(i+1,j,2)^2)+...
                 sqrt(vel_flux(i+1,j+1,1)^2+vel_flux(i+1,j+1,2)^2)+...
                 sqrt(vel_flux(i,j+1,1)^2+vel_flux(i,j+1,2)^2))/4;
    end
end

%       p -> contourplot

figure
contourf(Xc,Yc,P_ex); % contourlines for the pressure (exact solution)
axis equal		
% shading flat;
grid off		
set(gca, 'YTick', []);	
set(gca, 'XTick', []);
colorbar
title('Darcy Pressure:sol exact')

figure
contourf(Xc,Yc,Pc); % contourlines for the pressure (computed solution)
axis equal		
% shading flat;
grid off		
set(gca, 'YTick', []);	
set(gca, 'XTick', []);
colorbar
title('Darcy Pressure:computed sol')

%       u -> contourplot

figure
contourf(Xc,Yc,U_ex); % contourlines for the displacement (exact solution)
axis equal		
% shading flat;
grid off		
set(gca, 'YTick', []);	
set(gca, 'XTick', []);
colorbar
title('Displacement (Euclidean norm):sol exact')

figure
contourf(Xc,Yc,Uc); % contourlines for the displacement (computed solution)
axis equal		
% shading flat;
grid off		
set(gca, 'YTick', []);	
set(gca, 'XTick', []);
colorbar
title('Displacement (Euclidean norm):computed sol')

%       gamma -> contourplot

figure
contourf(Xc,Yc,G_ex); % contourlines for the rotation (exact solution)
axis equal		
% shading flat;
grid off		
set(gca, 'YTick', []);	
set(gca, 'XTick', []);
colorbar
title('Rotation: sol exact')

figure
contourf(Xc,Yc,Gc); % contourlines for the rotation (computed solution)
axis equal		
% shading flat;
grid off		
set(gca, 'YTick', []);	
set(gca, 'XTick', []);
colorbar
title('Rotation: computed sol')

%       z -> contourplot

figure
contourf(Xc,Yc,Z_ex); % contourlines for the velocity (exact solution)
axis equal		
% shading flat;
grid off		
set(gca, 'YTick', []);	
set(gca, 'XTick', []);
colorbar
title('Darcy velocity (Euclidean norm):sol exact')

figure
contourf(Xc,Yc,Zc); % contourlines for the velocity (computed solution)
axis equal		
% shading flat;
grid off		
set(gca, 'YTick', []);	
set(gca, 'XTick', []);
colorbar
title('Darcy velocity (Euclidean norm):computed sol')

%       sigma (1st row)-> contourplot

figure
contourf(Xc,Yc,S_ex1); % contourlines for the 1st row of the stress (exact solution)
axis equal		
% shading flat;
grid off		
set(gca, 'YTick', []);	
set(gca, 'XTick', []);
colorbar
title('Stress (1st row):sol exact')   
    
figure
contourf(Xc,Yc,Sc1); % contourlines for the 1st row of the stress (computed solution)
axis equal		
% shading flat;
grid off		
set(gca, 'YTick', []);	
set(gca, 'XTick', []);
colorbar
title('Stress (1st row):computed sol')

figure
contourf(Xc,Yc,S_ex2); % contourlines for the 2nd row of the stress (exact solution)
axis equal		
% shading flat;
grid off		
set(gca, 'YTick', []);	
set(gca, 'XTick', []);
colorbar
title('Stress (2nd row):sol exact')    
    
figure
contourf(Xc,Yc,Sc2); % contourlines for the 2nd row of the stress (computed sol)
axis equal		
% shading flat;
grid off		
set(gca, 'YTick', []);	
set(gca, 'XTick', []);
colorbar
title('Stress (2nd row):computed sol')
toc
end
function MS_MF_MFE(N)

global NN x y alpha lambda mu perm %c0

NN=N;       % Dimensión del problema discreto.

mesh=0;                 
[x,y]=init_mesh(mesh);  % coordenadas de los vértices de la malla.

% Parámetros ecuación Biot
alpha=0.93;
E=10^5;
nu=0.4;
lambda=E*nu/((1 + nu)*(1 - 2*nu));
mu=E/(2*(1 + nu));

% Storativity coefficient
% c0=0;

% Hydraulic conductivity: is inside the function 'kinv.m'
% K=[perm 0;0 perm];
perm=10^-7;

% initial time
t=0;
% Time step
delta_t=1e-03;
% Final time
Tf=1e-03;
% Tf=5*10^-3;
% Tf=1;

%  Matrices del sistema de Biot: A11, A12, A21 y A22
    % Asp y App las utilizaremos después
[A11,A12,A22,AspT,App,f_hat,q_hat]=build_matrices_and_vectors_Biot(delta_t);
A21=-A12';
Asp=AspT';

% Matriz del sistema reducido de Biot
Biot_matrix=[A11 A12;A21 A22];

%% Terms involving time

% Initial solution of the variables at t=0 
u=zeros(2*N*N,1);
p=zeros(N*N,1);

% for j=1:N
%     for i=1:N
%         ind2u=(i+(j-1)*N)*2;
%         ind1u=ind2u-1;
%         ind1p=ind2u/2;
%         
%         xx=(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1))/4;
%         yy=(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1))/4;
%         
%         u(ind1u)=sol_exactax(xx,yy,t,1);
%         u(ind2u)=sol_exactax(xx,yy,t,2);
%         p(ind1p)=sol_exactax(xx,yy,t,3);
%     end
% end

gamma=compute_gamma(u,p);
sigma=compute_tensors(u,p,gamma);

% % Initialize errors
% erroru_L2_inf=0;
% erroru_3_inf=erroru_L2_inf;
% errorg_L2_inf=erroru_L2_inf;
% error_sigma_sigmah_inf=erroru_L2_inf;
% errorp_L2_inf=erroru_L2_inf;
% errorp_3_inf=erroru_L2_inf;
% error_z_zh_inf=erroru_L2_inf;

% We solve Biot's system: Time loop

while t < Tf
    
    % Initial time terms affecting q
gTp=Asp*sigma + App*p;

q_term=q_hat + gTp;
    
    % Right-hand side of the Biot system
indep_term=[f_hat;q_term];
    
% Solution of the Biot system for the displacement and pressure vectors
sol_vec=Biot_matrix\indep_term; 
u=sol_vec(1:2*N*N);
p=sol_vec(2*N*N+1:2*N*N+N*N);

t=t+delta_t;

    % Now we compute the rest of the variables at the new time step t:
        % rotation 
gamma=compute_gamma(u,p);  % Computed solution for the rotation term
        % stress
sigma=compute_tensors(u,p,gamma);
%         % velocity
% [z,zx,zy]=compute_fluxes(p,t); 

%     % We reorder the variables to compute the errors and some contourplots
% gamma_n=reshape(gamma,N+1,N+1);
% sigma_n=build_sigma_n_2(sigma);
% vel_n=build_vel_n(z);

%     % L2 norms of the different variables
% % t
%         % Displacement
% [erroru_L2,erroru_3]=compute_error_displacements(u,t);
%         % Pressure
% [errorp_L2,errorp_3]=compute_error_pressures(p,t);
%         % Rotation
% errorg_L2=compute_error_rotations(gamma_n,t);
%         % Stress
% error_sigma_sigmah=compute_error_stress(sigma_n,nu,t);
%         % Velocity
% error_z_zh=compute_error_velocities(vel_n,t);
% 
%     % Infinity norm of the L2 errors of the variables
% erroru_L2_inf=max(erroru_L2_inf,erroru_L2);
% erroru_3_inf=max(erroru_3_inf,erroru_3);
% errorp_L2_inf=max(errorp_L2_inf,errorp_L2);
% errorp_3_inf=max(errorp_3_inf,errorp_3);
% errorg_L2_inf=max(errorg_L2_inf,errorg_L2);
% error_sigma_sigmah_inf=max(error_sigma_sigmah_inf,error_sigma_sigmah);
% error_z_zh_inf=max(error_z_zh_inf,error_z_zh);
end

% % We display the infinity norms of the errors of the variables at final time step
% erroru_L2_inf
% erroru_3_inf
% errorp_L2_inf
% errorp_3_inf
% errorg_L2_inf
% error_sigma_sigmah_inf
% error_z_zh_inf
% 
% % We display the L2 norms of the errors of the variables at final time step
% t
% erroru_L2
% erroru_3
% errorp_L2
% errorp_3
% errorg_L2
% error_sigma_sigmah
% error_z_zh

% disp(Biot_matrix)
% disp(indep_term)

% We compute the variable z of velocity at the final time step Tf(t):
% [z,zx,zy]=compute_fluxes(p,t); 

% Components of the displacement vector: u1 and u2
% u1=u(1:2:2*N*N);
% u2=u(2:2:2*N*N);

%% Plots of the Magnitude/Vectors of the variables

%% Plots for the pressure variable 

% We have to make N=10 for the plots

P=reshape(p,N,N);

%       P -> colormap
figure 
pcolor(P);
colorbar
title('Darcy Pressure, t=0.001')

%       P -> plot for different x-lines
Yc=zeros(1,N);

i=1;
for j=1:N
        yy=(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1))/4;
        Yc(i,j)=yy;
end

figure
plot(Yc,P(2,:),'o-blue');
hold on
plot(Yc,P(3,:),'d-red');
hold on
plot(Yc,P(4,:),'s-yellow');
hold on
plot(Yc,P(5,:),'|-magenta');
hold on
plot(Yc,P(6,:),'>-green');

legend('x=.15','x=.25','x=.35','x=.45','x=.55')

title('Cantilever bracket problem, t=0.005')
xlabel('y-coordinate')
ylabel('Pressure')

set(gca, 'YTick', [-1.6,-1.4,-1.2,-1,-0.8,-0.6,-0.4,-0.2,...
                    0,0.2,0.4,0.6,0.8,1,1.2,1.4]);	
set(gca, 'XTick', [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]);

end
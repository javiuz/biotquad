function MS_MF_MFE(N)

global NN x y lambda mu alpha c0
 
NN=N;       % Dimensión del problema discreto.

% Generación de la malla 
mesh=0;                 
[x,y]=init_mesh(mesh);  % coordenadas de los vértices de la malla.

% Parámetros ecuación Biot
alpha=1;
lambda=0;
mu=1;
% lambda=123;
% mu=79.3;

% Storativity coefficient
c0=0;

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

% Matriz de fuerza: A_sigma_p
A_sp=(A_sp_T)';

% Matriz de fuerza: pressure-pressure
A_pp=(c0+alpha^2/(lambda+mu))*(1/(N^2))*eye(N*N);

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
sigma=build_sigma_0_cartesian(t,nu);

%% Initialize errors
erroru_L2_inf=0;
erroru_3_inf=erroru_L2_inf;
errorg_L2_inf=erroru_L2_inf;
error_sigma_sigmah_inf=erroru_L2_inf;
errorp_L2_inf=erroru_L2_inf;
errorp_3_inf=erroru_L2_inf;
error_z_zh_inf=erroru_L2_inf;
   
%% We solve Biot's system: Time loop

while t < Tf

indep2=build_indep_f(t+delta_t);

indep5=delta_t*build_indep_q(t+delta_t) + A_sp*sigma + A_pp*p;

[indep1,indep4]=dir_bc_Pg(t+delta_t);

indep=[indep1;indep2;zeros((N+1)*(N+1),1);indep4;indep5];
   
Sol=M_big\indep;

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

end
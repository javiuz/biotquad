function MS_MF_MFE(N)

global NN x y

NN=N;       % Dimensi�n del problema discreto.

% Generaci�n de la malla 
mesh=0;                 
[x,y]=init_mesh(mesh);  % coordenadas de los v�rtices de la malla.

% Par�metros ecuaci�n Biot
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

%  Matriz de fuerza: stress-stress
A_ss=build_matrix_sig_sig(lambda,mu);

% Matriz de fuerza (transpuesta): stress-displacement
A_su_T=build_matrix_sig_u;

% Matriz de fuerza: A_sigma_u
A_su=(A_su_T)';

% Matriz de fuerza (transpuesta): stress-rotation
A_sg_T=build_matrix_sig_gamma;

% Matriz de fuerza: A_sigma_gamma
A_sg=(A_sg_T)';

% Matriz de fuerza (transpuesta): stress-pressure
A_sp_T=build_matrix_p_tau(alpha,lambda,mu);

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

% gamma=(A_sg*(A_ss\A_sg_T))\(-A_sg*(A_ss\A_su_T)*u -A_sg*(A_ss\A_sp_T)*p);
% sigma=A_ss\(-A_su_T*u -A_sg_T*gamma -A_sp_T*p);

%% We compute initial data for sigma, u and gamma from the elasticity problem
    
M_red_elas=[A_su*(A_ss\A_su_T)-...
            ((A_su*(A_ss\A_sg_T))*((A_sg*(A_ss\A_sg_T))\(A_sg*(A_ss\A_su_T))))];

indep2=build_indep_f(t);
        
indep_red_elas=indep2+(A_su-...
    ((A_su*(A_ss\A_sg_T))*((A_sg*(A_ss\A_sg_T))\(A_sg))))*(A_ss\(-A_sp_T*p));

u_elas=M_red_elas\indep_red_elas;

gamma_elas=(A_sg*(A_ss\A_sg_T))\(A_sg*(A_ss\(-A_sp_T*p))-...
            (A_sg*(A_ss\A_su_T))*u_elas);
        
sigma_elas=A_ss\(-A_sp_T*p -A_su_T*u_elas -A_sg_T*gamma_elas);

% M_elas=[A_ss A_su_T A_sg_T;...
%         -A_su zeros(2*N*N,2*N*N+(N+1)*(N+1));...
%         -A_sg zeros((N+1)*(N+1),2*N*N+(N+1)*(N+1))];
%    
% indep2=build_indep_f(t);
% 
% indep_elas=[-A_sp_T*p;indep2;zeros((N+1)*(N+1),1)];
% 
% Sol_elas=M_elas\indep_elas;
% 
% sigma_elas=Sol_elas(1:8*N*(N+1),1);
% u_elas=Sol_elas(8*N*(N+1)+1:8*N*(N+1)+2*N*N,1);
% gamma_elas=Sol_elas(8*N*(N+1)+2*N*N+1:end,1);

max(abs(sigma_elas-A_ss\(-A_sp_T*p -A_su_T*u_elas -A_sg_T*gamma_elas)))

indep5=delta_t*build_indep_q(t) + A_sp*sigma_elas + A_pp*p;

indep=[zeros(8*N*(N+1),1);indep2;zeros((N+1)*(N+1),1);...
       zeros(4*N*(N+1),1);indep5];
   
% indep_st=[zeros(8*N*(N+1),1);indep2;zeros((N+1)*(N+1),1);...
%        zeros(4*N*(N+1),1);delta_t*build_indep_q(t)];
   
Sol=M_big\indep;

% Sol_st=M_big_st\indep_st;

end
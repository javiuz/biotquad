global lambda mu alpha c0 

% Parámetros ecuación Biot
alpha=1;
lambda=1;
mu=1;

% Storativity coefficient
% c0=0;
c0=1e-05;

% Time step
t=0;
delta_t=1e-04;
Tf=1e-03;

%  Matriz de fuerza: stress-stress
A_ss=build_matrix_sig_sig;

% Matriz de fuerza (transpuesta): stress-displacement
A_su_T=build_matrix_sig_u;

% Matriz de fuerza: stress-displacement
A_su=(A_su_T)';

% Matriz de fuerza (transpuesta): stress-rotation
A_sg_T=build_matrix_sig_gamma;

% Matriz de fuerza: stress-rotation
A_sg=(A_sg_T)';

% Matriz de fuerza (transpuesta): stress-pressure
A_sp_T=build_matrix_sig_p;

% Matriz de fuerza: stress-pressure
A_sp=(A_sp_T)';

% Matriz de fuerza: velocity-velocity
A_zz=build_matrix_z_z;

% Matriz de fuerza (transpuesta): velocity-pressure
A_zp_T=build_matrix_z_p;

% Matriz de fuerza: velocity-pressure
A_zp=(A_zp_T)'*delta_t;

% Matriz de fuerza: pressure-pressure
A_pp=(c0+alpha^2/(lambda+mu))*(1/4)*eye(4);

% Mass matrix of the system 
MM=[A_ss A_su_T A_sg_T zeros(48,24) A_sp_T;...
    -A_su zeros(8,8+9+24+4);...
    -A_sg zeros(9,8+9+24+4);...
    zeros(24,48+8+9) A_zz A_zp_T;...
    A_sp zeros(4,8+9) -A_zp A_pp];

% Contribution of Dir. B.C of the displacement
Pgu=dir_bc_u;

% Initial values for the pressure and the stress
p=zeros(4,1);
sigma=build_initial_sigma;

% Flux source term q
q_term=(c0*delta_t/16)*[1;3;1;3];

while t <Tf

% Elasticity source term f
f=build_indep_f(t+delta_t);

% Contribution of Dir. B.C of the pressure
Pgp=dir_bc_p(t+delta_t);

% Vector q with all the contributions from previous time step
q=build_indep_q(q_term,sigma,p);

% Indep vector term
V_indep=[Pgu;f;zeros(9,1);Pgp;q];

% Solution of the system
Sol=MM\V_indep;

% New solution for the variables
sigma=Sol(1:48);
u=Sol(49:56);
gamma=Sol(57:65);
z=Sol(66:89);
p=Sol(90:93);

% Update t
t=t+delta_t;

end

sigma
u
gamma
z
p

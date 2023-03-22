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
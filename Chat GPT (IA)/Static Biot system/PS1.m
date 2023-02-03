% Define las constantes del problema
alpha = 1; % coeficiente de expansión térmica
lambda = 1; % módulo de Lamé
mu = 1; % módulo de Lamé
k = 1e-6; % conductividad térmica
c0 = 1e-5; % módulo de elasticidad volumétrica

% Define el dominio y la malla
x = [0, 1, 2, 3, 4, 5, 6, 7, 8];
y = [0, 1, 2, 3, 4, 5, 6, 7, 8];
[X, Y] = meshgrid(x, y);

% Inicializa las soluciones
u = zeros(size(X)); % deformación
p = zeros(size(X)); % presión

% Define las condiciones de contorno
u(1, :) = 1;
u(end, :) = 2;
p(:, 1) = 1;
p(:, end) = 2;

% Resuelve el sistema estático de Biot
for i = 2:length(x)-1
    for j = 2:length(y)-1
        % Calcula la deformación
        u(i, j) = alpha * p(i, j);
        
        % Calcula la presión
        p(i, j) = (c0 + lambda) * (u(i, j) - u(i-1, j)) + mu * (u(i, j) + u(i-1, j) - u(i, j-1) - u(i-1, j-1))/2/k;
    end
end

% Grafica las soluciones
surf(X, Y, u);
xlabel('x');
ylabel('y');
zlabel('Deformación');

figure;
surf(X, Y, p);
xlabel('x');
ylabel('y');
zlabel('Presión');

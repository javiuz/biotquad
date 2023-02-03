% Define las constantes del problema
alpha = 1; % coeficiente de expansión térmica
lambda = 1; % módulo de Lamé
mu = 1; % módulo de Lamé
k = 1e-6; % conductividad térmica
c0 = 1e-5; % módulo de elasticidad volumétrica
tau = 1e-4; % constante de tiempo

% Define el dominio y la malla
x = linspace(0, 1, 8);
y = linspace(0, 1, 8);
[X, Y] = meshgrid(x, y);

% Inicializa las soluciones
u = zeros(size(X)); % deformación
p = zeros(size(X)); % presión

% Define las condiciones de contorno
u(1, :) = 1;
u(end, :) = 2;
p(:, 1) = 1;
p(:, end) = 2;

% Define el tiempo final y la malla temporal
tiempo_final = 1e-3;
t = 0:tau:tiempo_final;

% Resuelve el sistema cuasi-estático de Biot en el tiempo
for n = 1:length(t)
    % Calcula la deformación
    for i = 2:length(x)-1
        for j = 2:length(y)-1
            u(i, j) = u(i, j) + tau * alpha * p(i, j);
        end
    end
    
    % Calcula la presión
    for i = 2:length(x)-1
        for j = 2:length(y)-1
            p(i, j) = p(i, j) + tau * (c0 + lambda) * (u(i, j) - u(i-1, j)) + tau * mu * (u(i, j) + u(i-1, j) - u(i, j-1) - u(i-1, j-1))/2/k;
        end
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

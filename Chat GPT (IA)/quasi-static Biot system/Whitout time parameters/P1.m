% Define las constantes del problema
alpha = 1; % coeficiente de expansión térmica
lambda = 1; % módulo de Lamé
mu = 1; % módulo de Lamé
k = 1e-6; % conductividad térmica

% Define el dominio y la malla
x = linspace(0, 1, 10);
y = linspace(0, 1, 10);
[X, Y] = meshgrid(x, y);

% Inicializa las soluciones
u = zeros(size(X)); % deformación
p = zeros(size(X)); % presión

% Define las condiciones de contorno
u(1, :) = 1;
u(end, :) = 2;
p(:, 1) = 1;
p(:, end) = 2;

% Resuelve el sistema cuasi-estático de Biot iterativamente
tol = 1e-6; % tolerancia
max_iter = 1000; % número máximo de iteraciones
err = inf; % error inicial
iter = 0;
while err > tol && iter < max_iter
    % Calcula la deformación
    for i = 2:length(x)-1
        for j = 2:length(y)-1
            u(i, j) = (u(i-1, j) + u(i+1, j) + u(i, j-1) + u(i, j+1))/4 + alpha * p(i, j);
        end
    end
    
    % Calcula la presión
    for i = 2:length(x)-1
        for j = 2:length(y)-1
            p(i, j) = (p(i-1, j) + p(i+1, j) + p(i, j-1) + p(i, j+1))/4 - (lambda + 2*mu) * (u(i, j) - u(i-1, j))/2/k;
        end
    end
    
    % Actualiza el error
    err = max(max(abs(u - u_prev)));
    iter = iter + 1;
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

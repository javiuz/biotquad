function [vel_flux]=compute_vel_fluxes(ux,uy)

global NN x y

nx=NN;
ny=NN;

% Componentes normales de las distintas fronteras.

normal_n=zeros(nx+1,ny+1,2);
normal_s=normal_n;
normal_e=normal_n;
normal_w=normal_n;

vel_flux=normal_n;

for i=1:nx+1
        j=1;
        normal_n(i,j,1)=(y(i,j+1)-y(i,j))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
        normal_n(i,j,2)=(x(i,j)-x(i,j+1))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
    for j=2:ny
        normal_n(i,j,1)=(y(i,j+1)-y(i,j))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
        normal_n(i,j,2)=(x(i,j)-x(i,j+1))/sqrt((y(i,j+1)-y(i,j))^2+(x(i,j)-x(i,j+1))^2);
        normal_s(i,j,1)=(y(i,j)-y(i,j-1))/sqrt((y(i,j)-y(i,j-1))^2+(x(i,j-1)-x(i,j))^2);
        normal_s(i,j,2)=(x(i,j-1)-x(i,j))/sqrt((y(i,j)-y(i,j-1))^2+(x(i,j-1)-x(i,j))^2);
    end
        j=ny+1;
        normal_s(i,j,1)=(y(i,j)-y(i,j-1))/sqrt((y(i,j)-y(i,j-1))^2+(x(i,j-1)-x(i,j))^2);
        normal_s(i,j,2)=(x(i,j-1)-x(i,j))/sqrt((y(i,j)-y(i,j-1))^2+(x(i,j-1)-x(i,j))^2);
end

for j=1:ny+1
        i=1;
        normal_e(i,j,1)=(y(i,j)-y(i+1,j))/sqrt((y(i,j)-y(i+1,j))^2+(x(i+1,j)-x(i,j))^2);
        normal_e(i,j,2)=(x(i+1,j)-x(i,j))/sqrt((y(i,j)-y(i+1,j))^2+(x(i+1,j)-x(i,j))^2);
    for i=2:nx
        normal_w(i,j,1)=(y(i-1,j)-y(i,j))/sqrt((y(i-1,j)-y(i,j))^2+(x(i,j)-x(i-1,j))^2);
        normal_w(i,j,2)=(x(i,j)-x(i-1,j))/sqrt((y(i-1,j)-y(i,j))^2+(x(i,j)-x(i-1,j))^2);
        normal_e(i,j,1)=(y(i,j)-y(i+1,j))/sqrt((y(i,j)-y(i+1,j))^2+(x(i+1,j)-x(i,j))^2);
        normal_e(i,j,2)=(x(i+1,j)-x(i,j))/sqrt((y(i,j)-y(i+1,j))^2+(x(i+1,j)-x(i,j))^2);
    end
        i=nx+1;
        normal_w(i,j,1)=(y(i-1,j)-y(i,j))/sqrt((y(i-1,j)-y(i,j))^2+(x(i,j)-x(i-1,j))^2);
        normal_w(i,j,2)=(x(i,j)-x(i-1,j))/sqrt((y(i-1,j)-y(i,j))^2+(x(i,j)-x(i-1,j))^2);
end

% Calculamos los vectores del flujo en cada uno de los nodos de la malla, 
% asociados según el tipo de celda considerada en cada caso.

for j=1:ny+1
    for i=1:nx+1
        % Flujo en la celda sureste, intervienen normal_s, normal_e
        denom=normal_s(i,j,1)*normal_e(i,j,2)-normal_e(i,j,1)*normal_s(i,j,2);
        if denom ~= 0
            alpha=normal_e(i,j,2)/denom;
            beta=-normal_s(i,j,2)/denom;
            gamma=-normal_e(i,j,1)/denom;
            delta=normal_s(i,j,1)/denom;
        vel_flux(i,j,1)=vel_flux(i,j,1)+alpha*ux(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
            beta*uy(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);
        vel_flux(i,j,2)=vel_flux(i,j,2)+gamma*ux(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2)+...
            delta*uy(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2);
        end
        % Flujo en la celda noreste, intervienen normal_e, normal_n
        denom=normal_e(i,j,1)*normal_n(i,j,2)-normal_n(i,j,1)*normal_e(i,j,2);
        if denom ~= 0
            alpha=normal_n(i,j,2)/denom;
            beta=-normal_e(i,j,2)/denom;
            gamma=-normal_n(i,j,1)/denom;
            delta=normal_e(i,j,1)/denom;
        vel_flux(i,j,1)=vel_flux(i,j,1)+alpha*uy(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
            beta*ux(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
        vel_flux(i,j,2)=vel_flux(i,j,2)+gamma*uy(i,j,1)/sqrt((x(i,j)-x(i+1,j))^2+(y(i,j)-y(i+1,j))^2)+...
            delta*ux(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
        end
        % Flujo en la celda noroeste, intervienen normal_n, normal_w
        denom=normal_n(i,j,1)*normal_w(i,j,2)-normal_w(i,j,1)*normal_n(i,j,2);
        if denom ~= 0
            alpha=normal_w(i,j,2)/denom;
            beta=-normal_n(i,j,2)/denom;
            gamma=-normal_w(i,j,1)/denom;
            delta=normal_n(i,j,1)/denom;
         vel_flux(i,j,1)=vel_flux(i,j,1)+alpha*ux(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
            beta*uy(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);
         vel_flux(i,j,2)=vel_flux(i,j,2)+gamma*ux(i,j,2)/sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2)+...
            delta*uy(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2);
        end
        % Flujo en la celda suroeste, intervienen normal_w, normal_s
        denom=normal_w(i,j,1)*normal_s(i,j,2)-normal_s(i,j,1)*normal_w(i,j,2);
        if denom ~= 0
            alpha=normal_s(i,j,2)/denom;
            beta=-normal_w(i,j,2)/denom;
            gamma=-normal_s(i,j,1)/denom;
            delta=normal_w(i,j,1)/denom;
        vel_flux(i,j,1)=vel_flux(i,j,1)+alpha*uy(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
            beta*ux(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);
        vel_flux(i,j,2)=vel_flux(i,j,2)+gamma*uy(i,j,2)/sqrt((x(i,j)-x(i-1,j))^2+(y(i,j)-y(i-1,j))^2)+...
            delta*ux(i,j,1)/sqrt((x(i,j)-x(i,j-1))^2+(y(i,j)-y(i,j-1))^2);
        end
    end
end

%%% Ponderación de los flujos en los distintos nodos %%%

% Las esquinas sólo tienen una celda asociada, no hace falta hacer
% ponderación.

% Los puntos de las fronteras tienen dos celdas asociadas, luego para 
% ponderar hay que dividir por dos.

vel_flux(2:nx,1,:)=vel_flux(2:nx,1,:)/2.; % Frontera Sur
vel_flux(1,2:ny,:)=vel_flux(1,2:ny,:)/2.; % Frontera Oeste
vel_flux(nx+1,2:ny,:)=vel_flux(nx+1,2:ny,:)/2.; % Frontera Este
vel_flux(2:nx,ny+1,:)=vel_flux(2:nx,ny+1,:)/2.; % Frontera Norte

% Los puntos interiores tienen cuatro celdas asociadas, luego para ponderar 
% hay que dividir por cuatro.

vel_flux(2:nx,2:ny,:)=vel_flux(2:nx,2:ny,:)/4.;

end
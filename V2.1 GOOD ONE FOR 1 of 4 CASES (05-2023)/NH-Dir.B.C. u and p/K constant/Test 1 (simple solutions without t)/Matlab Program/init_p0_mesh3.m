function p=init_p0_mesh3

global NN x y

N=NN;

% Initial solution of the pressure at t=0 
% p0 = @(x,y) y;
p=zeros(N*N,1);

for j=1:N
    for i=1:N
        ind1p=i+(j-1)*N;
        
        % Coordenadas de los vértices de la celda
        x1=x(i,j);
        y1=y(i,j);
        x2=x(i+1,j);
        y2=y(i+1,j);
        x3=x(i+1,j+1);
        y3=y(i+1,j+1);
        x4=x(i,j+1);
        y4=y(i,j+1);

        % Función a integrar en el elemento de referencia
        p0t = @(x,y) y1+(y2-y1).*x+(y4-y1).*y+(y3-y4-y2+y1).*x.*y;
            
        % Aproximación de la presión como la media de la presión en toda 
        % la celda (cuadrado unidad).
        p(ind1p) = integral2(p0t, 0, 1, 0, 1);
    end
end
return
end
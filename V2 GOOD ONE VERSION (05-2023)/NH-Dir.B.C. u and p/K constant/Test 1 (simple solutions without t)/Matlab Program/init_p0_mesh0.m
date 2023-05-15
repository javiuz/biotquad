function p=init_p0_mesh0

global NN x y

N=NN;

% Initial solution of the pressure at t=0 
p0 = @(x,y) y;
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
        
        xcord=[x1,x2,x3,x4];
        ycord=[y1,y2,y3,y4];
        
        % Límites de la celda rectangular
        xmin=min(xcord);
        ymin=min(ycord);
        xmax=max(xcord);
        ymax=max(ycord);
        
        % Área total de la celda rectangular
        Ac=area_cuadrilatero(x1,y1,x2,y2,x3,y3,x4,y4);
        
        % Aproximación de la presión como la media de la presión en toda 
        % la celda
        p(ind1p)=(1/Ac)*integral2(p0, xmin, xmax, ymin, ymax);
    end
end
return
end
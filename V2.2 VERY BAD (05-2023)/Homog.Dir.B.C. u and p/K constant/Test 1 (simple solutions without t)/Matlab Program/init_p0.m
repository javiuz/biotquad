function p=init_p0

global NN x y

N=NN;

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
        
        % Área total de la celda rectangular
        Ac=area_cuadrilatero(x1,y1,x2,y2,x3,y3,x4,y4);
        
        % Aproximación de la presión como la media de la presión en toda 
        % la celda
        p(ind1p)=(1/Ac)*(sol_exactax(x1,y1,0,3)+sol_exactax(x2,y2,0,3)+...
                         sol_exactax(x3,y3,0,3)+sol_exactax(x4,y4,0,3));
    end
end
return
end
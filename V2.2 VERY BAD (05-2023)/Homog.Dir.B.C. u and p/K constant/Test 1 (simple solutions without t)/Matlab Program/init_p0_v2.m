function p=init_p0_v2

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
        
        % Transformación de las coordenadas del punto medio del 
        % cuadrilátero a las del elemento de referencia
        xgc=x1 +(x2 - x1)*0.5 +(x4 - x1)*0.5 + (x3 - x4 - x2 + x1)*0.5^2;
        ygc=y1 +(y2 - y1)*0.5 +(y4 - y1)*0.5 + (y3 - y4 - y2 + y1)*0.5^2;
        
        %  Valor de p0(x,y) en el centro del cuadrilátero de referencia
        p(ind1p)=sol_exactax(xgc,ygc,0,3);
    end
end
return
end
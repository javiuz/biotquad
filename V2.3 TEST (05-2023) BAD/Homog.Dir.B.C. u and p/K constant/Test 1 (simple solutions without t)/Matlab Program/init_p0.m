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
       
        xc=(x1+x2+x3+x4)/4;
        yc=(y1+y2+y3+y4)/4;
        
        %  Interpolación del valor de P0(x,y) en el centro del cuadrilátero
        p(ind1p)=sol_exactax(xc,yc,0,3);
    end
end
return
end
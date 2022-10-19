function [erroru_L2,erroru_3]=compute_error_displacements(u,t)

global NN x y

N=NN;

erroru_L2=0;
erroru_3=0;

for j=1:N
    for i=1:N
        area=area_cuadrilatero(x(i,j),y(i,j),x(i+1,j),y(i+1,j),x(i+1,j+1),y(i+1,j+1),x(i,j+1),y(i,j+1));
        ind2=(i+(j-1)*N)*2;
        ind1=ind2-1;
        erroru_3=erroru_3+area*((sol_exacta(i,j,t,1)-u(ind1))^2+(sol_exacta(i,j,t,2)-u(ind2))^2);
        erroru_L2=erroru_L2+norma_2_gaussian(u(ind1:ind2),i,j,t);
    end
end

erroru_L2=sqrt(erroru_L2);
erroru_3=sqrt(erroru_3);

return
end
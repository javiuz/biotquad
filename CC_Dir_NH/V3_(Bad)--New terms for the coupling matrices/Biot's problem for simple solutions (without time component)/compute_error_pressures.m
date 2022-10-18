function [errorp_L2,errorp_3]=compute_error_pressures(p)

global NN x y

N=NN;

errorp_L2=0;
errorp_3=0;

for j=1:N
    for i=1:N
        area=area_cuadrilatero(x(i,j),y(i,j),x(i+1,j),y(i+1,j),x(i+1,j+1),y(i+1,j+1),x(i,j+1),y(i,j+1));
        errorp_3=errorp_3+area*(sol_exacta_p(i,j)-p(i+(j-1)*N))^2;
        errorp_L2=errorp_L2+norma_2_gaussian_p(p(i+(j-1)*N),i,j);
    end
end

errorp_L2=sqrt(errorp_L2);
errorp_3=sqrt(errorp_3);

return
end
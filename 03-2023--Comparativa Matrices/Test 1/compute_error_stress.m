function error_sigma_sigmah=compute_error_stress(sigma_n,nu,t)
global NN

nx=NN;
ny=NN;

error_sigma_sigmah=0;
for j=1:ny
    for i=1:nx
        intg=norma_2_gaussian_stress(i,j,sigma_n(i,j,:),nu,t);
        error_sigma_sigmah=error_sigma_sigmah+intg;
    end
end
error_sigma_sigmah=sqrt(error_sigma_sigmah);
return
end

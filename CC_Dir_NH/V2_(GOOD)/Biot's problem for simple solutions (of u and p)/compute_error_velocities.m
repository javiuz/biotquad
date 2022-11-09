function error_z_zh=compute_error_velocities(vel_n,t)
global NN

nx=NN;
ny=NN;

error_z_zh=0;
for j=1:ny
    for i=1:nx
        intg=norma_2_gaussian_vel(i,j,vel_n(i,j,:),t);
        error_z_zh=error_z_zh+intg;
    end
end
error_z_zh=sqrt(error_z_zh);
return
end
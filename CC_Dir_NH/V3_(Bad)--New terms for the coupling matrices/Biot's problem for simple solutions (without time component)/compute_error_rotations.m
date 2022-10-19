function errorg_L2=compute_error_rotations(g,t)

global NN

N=NN;

errorg_L2=0;

for j=1:N
    for i=1:N
        g1=g(i,j);
        g2=g(i+1,j);
        g3=g(i+1,j+1);
        g4=g(i,j+1);
        errorg_L2=errorg_L2+norma_2_gauss_lobatto_rot(g1,g2,g3,g4,i,j,t);
    end
end

errorg_L2=sqrt(2*errorg_L2);


return
end
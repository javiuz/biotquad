function App=build_matrix_p_p

global NN x y alpha lambda mu 

N=NN;

App=sparse(N*N,N*N);

coef_p=(alpha^2)/(16*mu*(lambda + mu));

% South-West corner node
i=1;
j=1;
ind1p=i+(j-1)*N;

x1=x(i,j);  
y1=y(i,j);            
x2=x(i+1,j);
y2=y(i+1,j);
x3=x(i+1,j+1);
y3=y(i+1,j+1);
x4=x(i,j+1);
y4=y(i,j+1);

Jer1=abs(x4*(y1 - y2) + x1*(y2 - y4) + x2*(-y1 + y4));
Jer2=abs(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3));
Jer3=abs(x4*(y2 - y3) + x2*(y3 - y4) + x3*(-y2 + y4));
Jer4=abs(x4*(y1 - y3) + x1*(y3 - y4) + x3*(-y1 + y4));

denom_p_1=coef_p/Jer1;
denom_p_2=coef_p/Jer2;
denom_p_3=coef_p/Jer3;
denom_p_4=coef_p/Jer4;

% Matriz local (A p, p)
k=term_App_matrix(x1,y1,x2,y2,x3,y3,x4,y4,denom_p_1,denom_p_2,denom_p_3,denom_p_4);

App(ind1p,ind1p)=k;

% South nodes (j=1)

for i=2:N
%     ind1p=i-1+(j-1)*N;
    ind2p=i+(j-1)*N;
    
%     x1=x(i-1,j);
%     y1=y(i-1,j);
    x2=x(i,j);
    y2=y(i,j);
    x3=x(i,j+1);
    y3=y(i,j+1);
%     x4=x(i-1,j+1);
%     y4=y(i-1,j+1);
    x5=x(i+1,j);
    y5=y(i+1,j);
    x6=x(i+1,j+1);
    y6=y(i+1,j+1);
    
    JE2r1=abs(x5*(-y2 + y3) + x3*(y2 - y5) + x2*(-y3 + y5));
    denom_p_1=coef_p/JE2r1;

    JE2r2=abs(x6*(y2 - y5) + x2*(y5 - y6) + x5*(-y2 + y6));
    denom_p_2=coef_p/JE2r2;
    JE2r3=abs(x6*(y3 - y5) + x3*(y5 - y6) + x5*(-y3 + y6));
    denom_p_3=coef_p/JE2r3;
    JE2r4=abs(x6*(-y2 + y3) + x3*(y2 - y6) + x2*(-y3 + y6));
    denom_p_4=coef_p/JE2r4;
    
%     k1=coef1_p;
    k2=term_App_matrix(x2,y2,x5,y5,x6,y6,x3,y3,denom_p_1,denom_p_2,denom_p_3,denom_p_4);
%     k=[0 0;0 k2];

    App(ind2p,ind2p)=k2;
end


for j=2:N

    % West nodes
    i=1;
    
%     ind1p=i+(j-2)*N;
    ind2p=i+(j-1)*N;
    
%     x1=x(i,j-1);
%     y1=y(i,j-1);
%     x2=x(i+1,j-1);
%     y2=y(i+1,j-1);
    x3=x(i+1,j);
    y3=y(i+1,j);
    x4=x(i,j);
    y4=y(i,j);
    x5=x(i+1,j+1);
    y5=y(i+1,j+1);
    x6=x(i,j+1);
    y6=y(i,j+1);

    JE2r1=abs(x6*(-y3 + y4) + x4*(y3 - y6) + x3*(-y4 + y6));
    denom_p_1=coef_p/JE2r1;
    
    JE2r2=abs(x5*(-y3 + y4) + x4*(y3 - y5) + x3*(-y4 + y5));
    denom_p_2=coef_p/JE2r2;
    JE2r3=abs(x6*(y3 - y5) + x3*(y5 - y6) + x5*(-y3 + y6));
    denom_p_3=coef_p/JE2r3;
    JE2r4=abs(x6*(y4 - y5) + x4*(y5 - y6) + x5*(-y4 + y6));
    denom_p_4=coef_p/JE2r4;
    
%     k1=coef1_p;
    k2=term_App_matrix(x4,y4,x3,y3,x5,y5,x6,y6,denom_p_1,denom_p_2,denom_p_3,denom_p_4);
%     k=[0 0;0 k2];
    
    App(ind2p,ind2p)=k2;
    
    % Central nodes 
         
    for i=2:N
%         ind1p=i-1+(j-2)*N;
%         ind2p=i+(j-2)*N;
%         ind3p=i-1+(j-1)*N;
        ind4p=i+(j-1)*N;
        
%     x1=x(i-1,j-1);
%     y1=y(i-1,j-1);
%     x2=x(i,j-1);
%     y2=y(i,j-1);
    x3=x(i,j);
    y3=y(i,j);
%     x4=x(i-1,j);
%     y4=y(i-1,j);
%     x5=x(i+1,j-1);
%     y5=y(i+1,j-1);
    x6=x(i+1,j);
    y6=y(i+1,j);
    x7=x(i+1,j+1);
    y7=y(i+1,j+1);
    x8=x(i,j+1);
    y8=y(i,j+1);
%     x9=x(i-1,j+1);
%     y9=y(i-1,j+1);
    
    JE3r1=abs(x8*(y3 - y6) + x3*(y6 - y8) + x6*(-y3 + y8));
    denom_p_1=coef_p/JE3r1;
    
    JE3r2=abs(x7*(y3 - y6) + x3*(y6 - y7) + x6*(-y3 + y7));
    denom_p_2=coef_p/JE3r2;
    JE3r3=abs(x8*(y6 - y7) + x6*(y7 - y8) + x7*(-y6 + y8));
    denom_p_3=coef_p/JE3r3;
    JE3r4=abs(x8*(y3 - y7) + x3*(y7 - y8) + x7*(-y3 + y8));
    denom_p_4=coef_p/JE3r4;
    
%     k1=coef1_p;
%     k2=coef2_p;
%     SIGUE EL ORDEN DE LAS PRESIONES EN EL SISTEMA: p4 va antes que p3
%     k4=coef4_p;
    k3=term_App_matrix(x3,y3,x6,y6,x7,y7,x8,y8,denom_p_1,denom_p_2,denom_p_3,denom_p_4);
%     k=[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 k3];
    
    App(ind4p,ind4p)=k3;
    end  
end
return
end
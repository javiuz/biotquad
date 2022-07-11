function k=kinv(r,c,i,j,k,l)

global x y NN

% Hydraulic conductivity
% K=[(xx+1)^2+yy^2 sin(xx*yy);sin(xx*yy) (xx+1)^2];

xx=(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1))/4;
yy=(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1))/4;
tensor=[(xx+1)^2+yy^2 sin(xx*yy);sin(xx*yy) (xx+1)^2];

% tensor_mean=tensor;

% if (i<=NN/2)
%     tensor=[1 0;0 1];
% else
%     tensor=1000*[1 0;0 1];
% end


% We compute the areas of triangles t1, t2 and t4 using the corresponding vertices
%coordinates
t1=area_triangulo(x(i,j+1),y(i,j+1),x(i,j),y(i,j),x(i+1,j),y(i+1,j));
t2=area_triangulo(x(i,j),y(i,j),x(i+1,j),y(i+1,j),x(i+1,j+1),y(i+1,j+1));
t4=area_triangulo(x(i+1,j+1),y(i+1,j+1),x(i,j+1),y(i,j+1),x(i,j),y(i,j));

if (k==i) && (l==j) % We are on the vertex (0,0) of the reference square
        xref=0.;
        yref=0.;
elseif (k==i+1) && (l==j)  % We are on the vertex (1,0) of the reference square
        xref=1.;
        yref=0.;
elseif (k==i+1) && (l==j+1)  % We are on the vertex (1,1) of the reference square
        xref=1.;
        yref=1.;
else % We are on the vertex (0,1) of the reference square
        xref=0.;
        yref=1.;
end 
df=[x(i+1,j)-x(i,j) x(i,j+1)-x(i,j);y(i+1,j)-y(i,j) y(i,j+1)-y(i,j)]+...
[(x(i+1,j+1)-x(i,j+1)-x(i+1,j)+x(i,j))*yref (x(i+1,j+1)-x(i,j+1)-x(i+1,j)+x(i,j))*xref;...
 (y(i+1,j+1)-y(i,j+1)-y(i+1,j)+y(i,j))*yref (y(i+1,j+1)-y(i,j+1)-y(i+1,j)+y(i,j))*xref];

jacob=2.*t1+2*(t2-t1)*xref+2*(t4-t1)*yref;
k_transf=jacob*inv(df)*tensor*(inv(df))';
kinv_matrix=inv(k_transf);
k=kinv_matrix(r,c);
return
end 
function sol=sol_exacta(i,j,~,comp)

global x y

xx=(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1))/4;
yy=(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1))/4;

% u=[u1;u2];

if comp==1
    sol=(1 - xx)*xx*yy;
else
    sol=xx*(1 - yy)*yy;
end

return
end
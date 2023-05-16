function sol=sol_exacta(i,j,tt,comp)

global x y

xx=(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1))/4;
yy=(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1))/4;

% u=[u1;u2];

if comp==1
    sol=exp(tt)*(xx^2 + xx^3*yy^4 + cos(1 - yy)*sin((1 - xx)*(1 - yy)));
else
    sol=exp(tt)*((1 - yy)^2 + (1 - xx)^4*(1 - yy)^3 + cos(xx*yy)*sin(xx));
end

return
end
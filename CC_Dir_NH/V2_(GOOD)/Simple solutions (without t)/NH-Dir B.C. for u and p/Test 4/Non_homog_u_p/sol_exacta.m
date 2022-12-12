function sol=sol_exacta(i,j,tt,comp)

global x y

xx=(x(i,j)+x(i+1,j)+x(i+1,j+1)+x(i,j+1))/4;
yy=(y(i,j)+y(i+1,j)+y(i+1,j+1)+y(i,j+1))/4;

% u=[u1;u2];

%         u1(xx,yy) = xx;
%         u2(xx,yy) = yy;

if comp==1
    sol=xx;
%     sol=0;
else
    sol=yy;
%     sol=0;
end

return
end
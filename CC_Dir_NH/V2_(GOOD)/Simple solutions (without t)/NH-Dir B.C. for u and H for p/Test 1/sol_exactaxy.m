function sol=sol_exactaxy(xx,yy,~)

%       u(xx,yy,tt) = {xx, 
%                    yy}

%         u1(xx,yy) = xx;
%         u2(xx,yy) = yy;

 sol=[xx^2,yy^2]';

return
end
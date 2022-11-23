function sol=sol_exactaxy(xx,yy,~)

%       u(xx,yy,tt) = {xx, 
%                    yy}

%         u1(xx,yy) = xx;
%         u2(xx,yy) = yy;

 sol=[xx*(1-xx)*yy*(1-yy),xx*(xx-1)*yy*(1-yy)]';

return
end
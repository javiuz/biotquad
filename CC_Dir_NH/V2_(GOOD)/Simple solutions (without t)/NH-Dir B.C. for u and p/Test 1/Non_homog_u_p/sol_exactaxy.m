function sol=sol_exactaxy(xx,yy,tt)

%       u(xx,yy,tt) = {xx, 
%                    yy}

%         u1(xx,yy) = xx;
%         u2(xx,yy) = yy;

 sol=[xx,yy]';
% sol=[tt*xx,tt*yy]';

return
end